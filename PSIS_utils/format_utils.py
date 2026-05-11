from config import *


def decimal_to_binary(n: int, bits: int = 8) -> list[int]:
    """
    将十进制转二进制(默认8位)
    :param n: 十进制整数
    :param bits: 要转换位数（默认8位
    :return: 二进制列表
    """
    n = n & ((1 << bits) - 1)  # 截断确保不超位数
    return [int(b) for b in format(n, f'0{bits}b')]


def binary_to_dict(binary: list[int]) -> dict:
    """
    二进制list转为多项式字典（0位均舍弃）
    :param binary: 二进制列表
    :return: 多项式 dict{}
    """
    reversed_binary = binary[::-1]
    return {i: coeff for i, coeff in enumerate(reversed_binary) if coeff != 0}


def binary_to_decimal(binary_list: list[int]) -> int:
    """
    二进制list转换为十进制数
    :param binary_list: 二进制list
    :return: 对应十进制
    """
    # 将列表元素转换为字符串并连接
    binary_str = ''.join(str(bit) for bit in binary_list)
    # 使用 int() 将二进制字符串转换为十进制数
    return int(binary_str, 2)


def poly_to_str(poly: dict) -> str:
    """
    美观展示，将多项式字典转为可读形式输出
    :param poly: 待展示多项式
    :return: 形式为：x^13 + x^10 + ...
    """
    terms = []
    for exp in sorted(poly.keys(), reverse=True):
        if exp == 0:
            terms.append("1")
        else:
            terms.append(f"x^{exp}" if exp > 1 else "x")
    return " + ".join(terms) if terms else "0"


def poly_to_bin(poly_dict: dict, bits: int) -> list:
    """
    将多项式转为指定位数的binary_list
    :param poly_dict: 多项式 dict{}
    :param bits: 需要转的二进制列表长度
    :return: 二进制列表
    """
    binary = [0] * bits  # Initialize a list of '0's with the specified length

    for exponent, coefficient in poly_dict.items():
        if coefficient == 1 and exponent < bits:  # Ensure valid exponent within bounds
            binary[bits - 1 - exponent] = 1  # Set the corresponding bit to '1'

    return binary


def poly_to_decimal(poly: dict) -> int:
    """
    将多项式dict转换为十进制表示
    :param poly: 待转换多项式 dict{}
    :return: 十进制整数
    """

    decimal_number = 0
    for exp, coeff in poly.items():
        if coeff == 1:
            decimal_number += 2 ** exp
    return decimal_number


def id_to_poly(id: int) -> dict:
    """
    输入十进制id，返回对应多项式字典
    :param id: 十进制id
    :return: id对应多项式 dict{}
    """
    id_bin = decimal_to_binary(id)
    return binary_to_dict(id_bin)


def reshape_pixels(pixels, rows, cols):
    total = rows * cols
    trimmed = pixels[:total]  # 截掉多余的
    reshaped = trimmed.reshape((rows, cols))  # 转成二维

    return reshaped


if __name__ == '__main__':
    import numpy as np

    pixels = np.arange(20)  # 长度为 20 的一维数组
    rows, cols = 3, 4  # 想要的二维形状

    print("==TEST==")
    print(reshape_pixels(pixels, rows, cols))