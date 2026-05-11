import cv2
import numpy as np
from config import *
from exception_handler import exception_handler
import random
from PSIS_utils.cal_utils import mul_poly, add_poly, poly_mod
from PSIS_utils.format_utils import *


def param_bin_to_decimal(binary_list: list[int], keys_list) -> list[int]:
    """
    传入格式：last pixel在前，高位在前，first pixel在后，低位在后
    [23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0]
    将恢复的binary_list传入，模板key_list传入，返回十进制恢复出的像素如[23,80,193]
    :param keys_list: 要恢复出的位
    :param binary_list:
    :return:
    """
    # 确保输入是有效的二进制列表
    if not all(bit in [0, 1] for bit in binary_list):
        raise ValueError("Input must be a list containing only 0s and 1s.")

    bin_len = len(binary_list)

    # TODO 未恢复位置为0、1、random？
    modified_set = MODIFIED_PIXEL_SET[MODIFIED_PIXEL_OPTION]

    # 初始化结果列表
    decimal_values = []

    # 按每8位分割列表
    for i in range(0, bin_len, 8):
        random_overall = random.randint(0, 1)
        # 对应i次项
        byte_list = binary_list[i:i + 8]
        for u in range(i, i + 8):
            if len(binary_list) - u - 1 not in keys_list:
                if modified_set == 'WHITE':
                    byte_list[u - i] = 1
                elif modified_set == 'BLACK':
                    byte_list[u - i] = 0
                elif modified_set == 'RANDOM':
                    byte_list[u - i] = random.randint(0, 1)
                else:
                    byte_list[u - i] = random_overall

        decimal_value = sum(bit * (2 ** (7 - idx)) for idx, bit in enumerate(byte_list))
        decimal_values.append(decimal_value)

    return list(decimal_values)


def remove_matching_terms(poly1: dict, poly2: dict):
    """
    从第一个多项式中删除与第二个多项式相同次数范围内的项。

    参数:
    poly1 (dict): 第一个多项式，以次数为键，系数为值。
    poly2 (dict): 第二个多项式，以次数为键，系数为值。

    返回:
    dict: 更新后的第一个多项式。
    """
    # 获取第一个多项式和第二个多项式的最高次数

    max_degree_poly1 = max(poly1.keys())
    max_degree_poly2 = max(poly2.keys())
    if max_degree_poly2 >= max_degree_poly1:
        return poly1
    # 计算需要删除的次数范围
    minus_degree = abs(max_degree_poly1 - max_degree_poly2)
    keys_to_remove = []
    for k in poly2.keys():
        for i in range(minus_degree + 1):
            keys_to_remove.append(k + i)

    keys_to_remove = set(keys_to_remove)

    for key in keys_to_remove:
        del poly1[key]

    return poly1


def get_pattern(m, bits_per_pixel=8):
    """
    获取8 * m的pattern (e(x))
    :return: dict{}
    """
    pattern = {}
    for i in range(bits_per_pixel * m):
        pattern[i] = 1
    return pattern


def permute_pixels(bin_list, permute_list=PERMUTE_LIST):
    to_deal_list = bin_list[::-1]
    # print("置乱前：", bin_list)
    n = len(bin_list)
    if len(permute_list) != n:
        raise ValueError("置乱模板出错！")

    # 先创建一个空列表来存放结果
    result = [None] * n
    for new_idx, old_idx in enumerate(permute_list):
        if not (0 <= old_idx < n):
            raise IndexError(f"模板出错：template 中的目标索引 {old_idx} 越界")
        result[new_idx] = to_deal_list[old_idx]

    result = result[::-1]
    return result


@exception_handler
def invert_permute(bin_list, permute_list=PERMUTE_LIST):
    """
    计算给定置乱模板 template 的逆模板 inv，
    使得 inv[ template[i] ] = i。
    :param bin_list: 原始bin_list，格式为[23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0]
    :param permute_list: 置乱模板，config传入
    :return:
    """
    n = len(permute_list)
    inv = [-1] * n

    for i, j in enumerate(permute_list):
        inv[j] = i

    result = permute_pixels(bin_list, inv)

    return result


def scramble_image_with_seed(image: np.ndarray, seed: int) -> np.ndarray:
    """
    使用给定的随机种子对灰度图像进行像素级置乱。
    :param image: 输入的灰度图像
    :param seed: 用于置乱的随机种子
    :return: 置乱后的图像
    """

    h, w = image.shape
    flat_image = image.flatten()
    rng = np.random.default_rng(seed)
    perm = rng.permutation(h * w)
    scrambled = flat_image[perm]
    return scrambled.reshape((h, w))


def unscramble_image_with_seed(scrambled_image: np.ndarray, seed: int) -> np.ndarray:
    """
    使用相同的随机种子对置乱后的灰度图像进行逆置乱
    :param scrambled_image: 置乱后的灰度图像
    :param seed: 与置乱时相同的随机种子
    :return: 恢复后的原始图像
    """
    h, w = scrambled_image.shape
    flat_scrambled = scrambled_image.flatten()
    rng = np.random.default_rng(seed)
    perm = rng.permutation(h * w)
    inverse_perm = np.argsort(perm)
    restored = flat_scrambled[inverse_perm]
    return restored.reshape((h, w))


def arnold_scramble(matrix, iterations=1, a=1, b=1):
    """
    对二维矩阵进行 Arnold 置乱。

    参数：
        matrix (np.ndarray): 输入的二维矩阵，必须是方阵。
        iterations (int): 置乱的迭代次数。
        a (int): 置乱参数 a，默认值为 1。
        b (int): 置乱参数 b，默认值为 1。

    返回：
        np.ndarray: 置乱后的矩阵。
    """
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Arnold 置乱要求输入矩阵为方阵。")

    N = matrix.shape[0]
    result = np.copy(matrix)

    for _ in range(iterations):
        temp = np.zeros_like(result)
        for x in range(N):
            for y in range(N):
                new_x = (a * x + y) % N
                new_y = (b * x + (a * b + 1) * y) % N
                temp[new_x, new_y] = result[x, y]
        result = temp

    return list(result)


def arnold_unscramble(matrix, iterations=1, a=1, b=1):
    """
    对二维矩阵进行 Arnold 逆置乱。

    参数：
        matrix (np.ndarray): 输入的二维矩阵，必须是方阵。
        iterations (int): 逆置乱的迭代次数。
        a (int): 置乱参数 a，默认值为 1。
        b (int): 置乱参数 b，默认值为 1。

    返回：
        np.ndarray: 逆置乱后的矩阵。
    """
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Arnold 逆置乱要求输入矩阵为方阵。")

    N = matrix.shape[0]
    result = np.copy(matrix)

    for _ in range(iterations):
        temp = np.zeros_like(result)
        for x in range(N):
            for y in range(N):
                new_x = ((a * b + 1) * x - b * y) % N
                new_y = (-a * x + y) % N
                temp[new_x, new_y] = result[x, y]
        result = temp

    return result


def open_image(file_path):
    """
    打开指定路径文件，然后返回像素值二维矩阵
    :param file_path:
    :return:
    """
    # TODO 这里仅支持灰度图像，多通道图像也许未来也可以支持
    image = cv2.imread(file_path, cv2.IMREAD_GRAYSCALE)
    return image.tolist()


def save_share(file_name, pixels, add_p=''):
    array = np.array(pixels, dtype=np.uint8)
    cv2.imwrite(file_name, array)
    print(f"SAVE IMAGE:{add_p} {file_name} is saved.")


def main_poly_cal(params: dict, id: dict, Gx: dict):
    # h = get_highest_degree(params)
    _base_ = {0: 0}
    for key in params.keys():
        _bs_ = {0: 1}
        for i in range(int(key)):
            _bs_ = mul_poly(id, _bs_)
        _bs_ = mul_poly(params[key], _bs_)
        _base_ = add_poly(_bs_, _base_)

    _base_ = poly_mod(_base_, Gx)
    return _base_
