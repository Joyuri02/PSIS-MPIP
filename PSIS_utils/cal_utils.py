from collections import defaultdict
from itertools import combinations
from sympy import symbols, Poly, GF
from config import *


def add_poly(poly1: dict, poly2: dict) -> dict:
    """
    计算 poly1 + poly2
    :param poly1:
    :param poly2:
    :return: 多项式dict{}
    """
    result = poly1.copy()  # 先复制 poly1 的项
    for exp, coeff in poly2.items():
        if exp in result:
            result[exp] += coeff  # 同类项相加
        else:
            result[exp] = coeff
    # 可选：将系数为 0 的项删除
    result = {exp: coeff for exp, coeff in result.items() if coeff != 0}
    return result


def sub_poly(poly1: dict, poly2: dict) -> dict:
    """
    计算 poly1 - poly2
    :param poly1: 被减(减号前)
    :param poly2: 减项(减号后)
    :return: 多项式dict{}
    """
    result = poly1.copy()
    for exp, coeff in poly2.items():
        if exp in result:
            result[exp] -= coeff
        else:
            result[exp] = -coeff
    # 移除系数为零的项
    result = {exp: coeff for exp, coeff in result.items() if coeff != 0}
    return result


def mul_poly(poly1: dict, poly2: dict) -> dict:
    """
    计算 poly1 * poly2
    :param poly1:
    :param poly2:
    :return: 多项式dict{}
    """
    result = defaultdict(int)
    for exp1, coeff1 in poly1.items():
        for exp2, coeff2 in poly2.items():
            result[exp1 + exp2] += coeff1 * coeff2
    return {k: v for k, v in result.items() if v != 0}


def gf2_inverse(poly_dict: dict, modulus_poly: dict) -> dict:
    """
    计算poly_dict 在 mod modulus_poly下的逆元
    :param poly_dict:
    :param modulus_poly: 要模的多项式 dict{}
    :return: 逆元多项式 dict{}
    """
    x = symbols('x')
    # 将字典表示的多项式转换为 Poly 对象
    poly = Poly(sum(coeff * x ** exp for exp, coeff in poly_dict.items()),
                x, domain=GF(2, symmetric=False))

    # 确保模多项式也在 GF(2, symmetric=False) 上
    modulus_poly = Poly(modulus_poly.as_expr(), x, domain=GF(2, symmetric=False))

    # 扩展欧几里得算法求逆：
    # 返回值顺序为 (s, t, h) 使得 s*poly + t*modulus_poly = h
    inverse_poly, _, gcd = poly.gcdex(modulus_poly)

    if gcd.as_expr() != 1:
        return None  # 多项式不可逆

    # 对于一元多项式，指数是单元素元组，转换为字典形式
    return {exp[0]: coeff for exp, coeff in inverse_poly.terms() if coeff != 0}


def generate_gf2n_polynomials(n: int) -> list:
    """
    生成 GF(2^n) 中的所有多项式，每个多项式以字典形式表示。
    :param n: 多项式的最高次幂
    :return: 包含所有多项式的列表，每个多项式以字典形式表示。
    """
    poly_list = []
    # 遍历可能的项数（从1到n）
    for num_terms in range(1, n + 1):
        # 在 [0, n-1] 范围内选择 num_terms 个不同的指数组合
        for exponents in combinations(range(n), num_terms):
            # 创建一个多项式字典，键为指数，值为1（GF(2)中的系数）
            poly = {exp: 1 for exp in exponents}
            poly_list.append(poly)
    return poly_list


def get_highest_degree(poly_dict: dict) -> int:
    """
    获取多项式的最高次数
    :param poly_dict: 待求多项式
    :return: 最高次数n
    """

    if not poly_dict:
        return -1  # 空字典表示零多项式，次数定义为-1
    # TODO 这里也许会有特殊情况？
    return max(poly_dict.keys())


def find_inverses(poly_list: list[dict], mod_poly: dict) -> list[dict]:
    """
    计算poly_list中的各个poly在mod mod_poly下的逆元
    :param poly_list: 待求逆的多项式列表 [{},{},...]
    :param mod_poly: 模多项式 dict{}
    :return: 逆元多项式列表 [{}, {}, {}]
    """
    x = symbols('x')

    # 将字典转换为SymPy多项式
    mod_poly = sum(coeff * x ** exp for exp, coeff in mod_poly.items())
    modulus_poly = Poly(mod_poly, x, domain=GF(2, symmetric=False))

    # 验证是否为不可约多项式（可选）
    # if not modulus_poly.is_irreducible:
    #     raise ValueError("模多项式必须是GF(2)上的不可约多项式")

    return [gf2_inverse(p, modulus_poly) for p in poly_list]


def poly_mod_coeff(poly: dict, mod: int = 2) -> dict:
    """
    对多项式系数取模（处理负系数）
    :param poly: 原始多项式字典 {指数: 系数}
    :param mod: 模数（如2, 7, 256等）
    :return: 系数取模后的字典
    """
    result = {}
    for exp, coeff in poly.items():
        # 处理负数：例如-5 mod 3 = 1
        adjusted = coeff % mod
        if adjusted != 0:  # 忽略系数为0的项
            result[exp] = adjusted
    return result


def poly_mod(poly: dict, generator: dict) -> dict:
    """
    计算 GF(2)下 poly % generator
    计算多项式模生成多项式的余式
    :param poly: 被除多项式 dict{}
    :param generator: 生成多项式 dict{}
    :return: 余式 dict{}
    """
    dividend = poly.copy()
    gen_degree = max(generator.keys())
    while True:
        # 获取当前被除式的最高次项
        current_degree = max(dividend.keys()) if dividend else -1
        if current_degree < gen_degree:
            break
        # 计算商项：系数 = 当前最高项系数 / 生成式最高项系数
        lead_coeff = dividend[current_degree]
        # 生成式乘商项并做减法（模2时等效异或）
        term = {exp + (current_degree - gen_degree): coeff * lead_coeff
                for exp, coeff in generator.items()}
        # 合并到被除式（若系数模运算需在此处处理）
        for exp, coeff in term.items():
            if exp in dividend:
                dividend[exp] -= coeff
                if dividend[exp] == 0:
                    del dividend[exp]
            else:
                dividend[exp] = -coeff
    # return dividend
    return poly_mod_coeff(dividend, 2)


def fix_to_ceil(n: int, k: int) -> int:
    padded = ((n + k - 1) // k) * k
    return padded - n


if __name__ == '__main__':
    n = 65536
    k = 6
    print("==TEST==")
    print(f"fix_to_ceil({n},{k}) = {fix_to_ceil(n, k)}")