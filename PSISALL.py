import time

from PSIS_utils.format_utils import *
from PSIS_utils.cal_utils import *
from PSIS_utils.PSISutils import *
from image_share import *
# from config import *
from config_poly_group import ALL_GROUP


def get_inverses(polys):
    """
    :param polys: list[dict{}]
    :return: level_polys: list[], invs: list[dict{}]
    """
    pmpt = "全量"
    if MODE_SMALL_TEST:
        pmpt = "测试"
    print(f"======MODE: {pmpt}======")
    levels = len(polys)
    level_polys = []
    ids = get_highest_degree(polys[0])
    # print("ids =", ids)
    # 首先要根据级数和ID数创建二维数组 行是ID数+1，列是级数+1
    invs = [[{} for _ in range(levels + 1)] for _ in range(2 ** ids)]
    # TODO 第一个一定是最低次的
    min_deg = get_highest_degree(polys[0])
    ids_list = generate_gf2n_polynomials(min_deg)
    # print(ids_list)
    # 获取G(x)
    # 首先polys的长度应该是多项式个数，如果是三个，应该有[0]mod g3, [1]mod g2xg3, [2]mod Gx
    # 如果是四个，应该是[0]mod g4, [1]mod g4xg3, [2]mod g4xg3xg2, [3]mod Gx
    tmp_poly = {0: 1}
    level_polys.append(tmp_poly)
    for i in range(levels - 1, -1, -1):
        # 从gn到g1
        tmp_poly = mul_poly(tmp_poly, polys[i])
        tmp_poly = poly_mod_coeff(tmp_poly, 2)
        level_polys.append(tmp_poly)
        # 所有id在mod(tmp_poly)下的逆元
        tmp_invv = find_inverses(ids_list, tmp_poly)
        for id in range(1, len(ids_list) + 1):
            invs[id][levels - i] = tmp_invv[id - 1]

    return level_polys, invs


def lagarange_interpolation(pairs, mod_poly, invs, lv, k=K):
    if k == 2:
        return lagarange_interpolation_2(pairs, mod_poly, invs, lv)
    if k == 3:
        return lagarange_interpolation_3(pairs, mod_poly, invs, lv)


def lagarange_interpolation_2(pairs, mod_poly, invs, lv):
    """
    实现传入 2 个id值+多项式，返回一个LI多项式（这里多项式的形式应该是十进制了）
    :return:
    """
    # 这里每个pair考虑用元组存储，pair[0]是id值，pair[1]是多项式（share），piars整体是list
    # TODO 扩展k个
    ids = [id_to_poly(pair[0]) for pair in pairs]
    # 有 k 个id

    a0 = {0: 0}
    a1 = {0: 0}
    # all_poly = [{} * k]  # k 个项
    base_x = {1: 1}
    for i, id in enumerate(ids):
        # 处理id这个项
        fi = pairs[i][1]
        p0 = fi
        p1 = fi
        for j, jid in enumerate(ids):
            if i == j:
                continue
            p1 = mul_poly(p1, jid)
            imj = sub_poly(id, jid)
            imj = poly_mod(imj, mod_poly)  # 先mod，再查逆元
            imj_idx = poly_to_decimal(imj)
            imj_inv = invs[imj_idx][lv]  # 对应级别下的逆元
            p0 = mul_poly(imj_inv, p0)
            p1 = mul_poly(imj_inv, p1)

        a0 = add_poly(a0, p0)
        a0 = poly_mod(a0, mod_poly)
        a1 = add_poly(a1, p1)
        a1 = poly_mod(a1, mod_poly)
    param0 = poly_to_decimal(a0)
    param1 = poly_to_decimal(a1)

    return {1: param0, 0: param1}


def lagarange_interpolation_3(pairs, mod_poly, invs, lv):
    """
    实现传入 3 个id值+多项式，返回一个LI多项式
    :param pairs:
    :param mod_poly:
    :param invs:
    :param lv:
    :return:
    """
    ids = [id_to_poly(pair[0]) for pair in pairs]

    a0 = {0: 0}
    a1 = {0: 0}
    a2 = {0: 0}
    base_x = {1: 1}
    for i, id in enumerate(ids):
        # 处理id这个项
        fi = pairs[i][1]
        p0 = {0: 1}
        p1 = {0: 0}
        p2 = {0: 1}
        all_inv = {0: 1}

        for j, jid in enumerate(ids):
            if i == j:
                continue
            imj = sub_poly(id, jid)
            imj = poly_mod(imj, mod_poly)  # 先mod，再查逆元
            imj_idx = poly_to_decimal(imj)
            imj_inv = invs[imj_idx][lv]  # 对应级别下的逆元
            all_inv = mul_poly(all_inv, imj_inv)

            p0 = mul_poly(p0, jid)
            p1 = sub_poly(p1, jid)

        all_inv = poly_mod(all_inv, mod_poly)
        all_inv = mul_poly(fi, all_inv)
        all_inv = poly_mod(all_inv, mod_poly)

        p0 = mul_poly(p0, all_inv)
        p1 = mul_poly(p1, all_inv)
        p2 = mul_poly(p2, all_inv)

        p0 = poly_mod(p0, mod_poly)
        p1 = poly_mod(p1, mod_poly)
        p2 = poly_mod(p2, mod_poly)
        # print(a0, p0)
        a0 = add_poly(a0, p0)
        a1 = add_poly(a1, p1)
        a2 = add_poly(a2, p2)
    a0 = poly_mod(a0, mod_poly)
    a1 = poly_mod(a1, mod_poly)
    a2 = poly_mod(a2, mod_poly)

    param0 = poly_to_decimal(a0)
    param1 = poly_to_decimal(a1)
    param2 = poly_to_decimal(a2)

    return {0: param0, 1: param1, 2: param2}


def lagrange_interpolation_k(pairs, mod_poly, invs, lv):
    """
    未使用的方法
    """
    ids = [id_to_poly(pair[0]) for pair in pairs]  # 3 个 ID 多项式
    values = [pair[1] for pair in pairs]  # 3 个 share 多项式
    result_poly = {0: 0}

    for i in range(3):
        xi = ids[i]
        yi = values[i]
        numerator = {0: 1}
        denominator = {0: 1}

        for j in range(3):
            if i == j:
                continue
            xj = ids[j]

            # 分子累乘 (x - xj)
            numerator = mul_poly(numerator, xj)  # x - xj：用 jid 表示 xj，多项式形式
            # 分母 xi - xj，模后再取逆
            diff = sub_poly(xi, xj)
            diff = poly_mod(diff, mod_poly)
            idx = poly_to_decimal(diff)
            inv = invs[idx][lv]
            denominator = mul_poly(denominator, inv)

        # lj(x) = numerator * denominator
        lj = mul_poly(numerator, denominator)
        lj = poly_mod(lj, mod_poly)

        # term = y_i * lj(x)
        term = mul_poly(yi, lj)
        term = poly_mod(term, mod_poly)

        result_poly = add_poly(result_poly, term)
        result_poly = poly_mod(result_poly, mod_poly)

    return result_poly  # 一个完整的 dict，形如 {2: xx, 1: yy, 0: zz}


def get_params(r_list, poly, m=M):
    """
    希望实现输入dict{}LI插值结果，提取各个参数，按顺序返回各个像素
    r_list : 在这个范围内的保留，为揭示出来的位
    :return:
    """
    len_poly = len(poly)  # 应该是k
    # 一共len_poly * m个像素
    # TODO 这里对于modified bits的处理有待商榷
    pixels = []
    for i in range(len_poly):
        param = poly.get(i)
        if param == 0:
            pixels.extend([0] * m)
            continue
        param_bin = decimal_to_binary(param, m * 8)
        # print(param_bin)
        param_list = param_bin_to_decimal(param_bin, r_list)
        # print("param_list =", param_list)
        pixels.extend(param_list)
        # print("pixels =", pixels)
    return pixels


def config_reader(poly_name, levels):
    config_group = ALL_GROUP.get(levels)
    cfg = config_group.get(poly_name)
    g1 = cfg.get("G1")
    g2 = cfg.get("G2")
    g3 = cfg.get("G3")
    g4 = cfg.get("G4")
    g5 = cfg.get("G5")
    g6 = cfg.get("G6")
    g7 = cfg.get("G7")
    m = cfg.get("M")
    polys = [g1, g2, g3, g4, g5, g6, g7]
    polys = polys[:levels]
    permute_list = cfg.get("permute_list")

    return polys, m, permute_list


def PSIS_XY(polys, m, permuted, permute_list, levels, recovered_path):
    all_times = []
    # # 加密
    # g1 = CONFIG_G1
    # g2 = CONFIG_G2
    # g3 = CONFIG_G3
    # g4 = CONFIG_G4
    # polys = [g1, g2, g3, g4]
    # polys = polys[:LEVELS]
    # # if len(polys) != LEVELS:
    # #     raise ValueError('传入多项式个数错误！len(polys) != LEVELS')

    # 二级恢复
    level_polys, invs = get_inverses(polys)  # 获取了逆元列表
    print(invs)
    # for ik, iv in enumerate(invs):
    #     print(ik,":",[poly_to_str(v) for v in iv])
    # for l in level_polys:
    #     print(poly_to_str(l))
    # print(level_polys)  # level_polys的长度是LEVELS，level_polys[1]是g4

    # print("polys = ", polys)
    # print("permute_list =", permute_list)
    from pathlib import Path

    folder_path = Path(TEST_IMG_PATH)  # 例如 Path(r"E:\images")
    image_files = list(folder_path.glob("*.png")) + \
                  list(folder_path.glob("*.pgm")) + \
                  list(folder_path.glob("*.bmp")) + \
                  list(folder_path.glob("*.tif"))

    image_files = [str(f) for f in image_files]

    print(image_files)
    for img_path in image_files:  # 对每张secret
        img_info = img_path.split("\\")[-1].split('.')
        img_name = img_info[0]
        img_suffix = img_info[-1]
        print(f"=================================================Processing Image: {img_name}.{img_suffix}")
        share_name_p = f'p{img_name}_share'
        res_name_p = f'p{img_name}_recovered'

        image_pixels = open_image(img_path)
        rows = len(image_pixels)
        cols = len(image_pixels[0])
        id1 = 2
        id2 = 3
        id3 = 4
        ids = [id1, id2]
        if K == 3:
            ids = [id1, id2, id3]

        start = time.time()
        small_start = time.time()
        # 原始秘密图像置乱
        scrambled = None
        image_pixels_np = np.array(image_pixels)
        if SCRAMBLED_OPTION == 'RANDOM':
            scrambled = scramble_image_with_seed(image_pixels_np, SEED)
        if SCRAMBLED_OPTION == 'ARNOLD':
            scrambled = arnold_scramble(image_pixels_np, ARNOLD_ITERATIONS)
        # end = time.time()
        # print(f"置乱耗时：{end - small_start}")
        # small_start = time.time()
        # print(level_polys)
        # 获取share1和share2（这里没有获取子share）
        # shares_pixels_dict = image_to_share(scrambled, ids, level_polys[levels], m=m,
        #                                     permuted=permuted, permute_list=permute_list)

        shares_pixels_dict = image_to_share_new(scrambled, ids, level_polys[levels], m=m,
                                                permuted=permuted, permute_list=permute_list)
        reshape_shares(shares_pixels_dict, ids, 512, 256)
        # end = time.time()
        # print(f"生成share耗时：{end - small_start}")
        # small_start = time.time()
        """
        当K=3时，这里要换成171
        """

        # print("==============shares============")
        # print(shares_pixels_dict)
        # TODO shares整理
        # reshape_shares(shares_pixels_dict, ids, rows, m * (cols // (K * m)))

        shares = []
        for id in ids:
            share_name = share_name_p + f'{id}'  # pxxx_share2、pxxx_share3
            share_save_path = SHARE_PATH + share_name + SHARE_SUFFIX
            share_pixels = shares_pixels_dict[id]
            shares.append(share_name)
            save_share(share_save_path, share_pixels)

        # shares的长度（share的数量）应与N相等
        # shares = [share1_name, share2_name]
        # 生成sub share
        for share in shares:
            share_path = SHARE_PATH + share + SHARE_SUFFIX
            pixels_share = open_image(share_path)
            # sub_share的数量是LEVELS - 1
            for lv in range(1, levels):
                mod_poly = level_polys[lv]
                pixels_share_sub = get_sub_share(pixels_share, mod_poly, m=m, reshape_to=256)
                share_file_path = SHARE_PATH + share + f'_{lv}' + SHARE_SUFFIX
                save_share(share_file_path, pixels_share_sub)
        # end = time.time()
        # print(f"生成subshares耗时：{end - small_start}")

        # 解密
        shares_s = []
        for lv in range(1, levels):
            # start = time.time()
            small_start = time.time()
            if lv == levels:
                share1_path = SHARE_PATH + share_name_p + f'{id1}' + SHARE_SUFFIX
                share2_path = SHARE_PATH + share_name_p + f'{id2}' + SHARE_SUFFIX
                share3_path = SHARE_PATH + share_name_p + f'{id3}' + SHARE_SUFFIX
                res_name = res_name_p
            else:
                share1_path = SHARE_PATH + share_name_p + f'{id1}_{lv}' + SHARE_SUFFIX
                share2_path = SHARE_PATH + share_name_p + f'{id2}_{lv}' + SHARE_SUFFIX
                share3_path = SHARE_PATH + share_name_p + f'{id3}_{lv}' + SHARE_SUFFIX
                res_name = res_name_p + f'_{lv}'

            pixels_share1 = open_image(share1_path)
            # print(share2_path)
            pixels_share2 = open_image(share2_path)
            # print("===========sub1========")
            # print(pixels_share1)
            pss = [pixels_share1, pixels_share2]
            if K == 3:
                pixels_share3 = open_image(share3_path)
                pss = [pixels_share1, pixels_share2, pixels_share3]
            mod_poly = level_polys[lv]
            # 这里不包括完全恢复
            rr = remove_matching_terms(get_pattern(m), mod_poly)
            print("RECOVER_LEVEL =", lv)
            print("unmodified bits =", list(rr.keys()))
            # rec = share_recover_2(pss, ids, mod_poly, invs, rr.keys(), lv,
            #                       m=m, permuted=permuted, permute_list=permute_list)

            rec = share_recover_3(pss, ids, mod_poly, invs, rr.keys(), lv,
                                  m=m, permuted=permuted, permute_list=permute_list)
            # print(rec)
            rec = pixels_1D_to_2D(rec, 512, 512)

            rec_np = np.array(rec)

            # 解密阶段
            rec_uncrambled = None
            if SCRAMBLED_OPTION == 'RANDOM':
                rec_uncrambled = unscramble_image_with_seed(rec_np, SEED)
            if SCRAMBLED_OPTION == 'ARNOLD':
                rec_uncrambled = arnold_unscramble(rec_np, ARNOLD_ITERATIONS)

            res_path = recovered_path + res_name + RECOVERD_SUFFIX
            save_share(res_path, rec_uncrambled)
            # end = time.time()
            # print(f"恢复耗时 {end - small_start}")
        end = time.time()
        print(f"总耗时 {end - start}")
        all_times.append(end - start)
    average = sum(all_times) / len(all_times)
    print("平均耗时:", average)


if __name__ == '__main__':
    permuted = False
    # fname = "Fy"
    fnames = [
        "F10"
    ]
    levels = 2
    for fname in fnames:
        polys, m, permute_list = config_reader(fname, levels)
        recovered_path = RECOVERED_PREFIX + "/" + fname + "/"
        if not permuted:
            recovered_path = RECOVERED_PREFIX + "/" + fname + "_w/"
        # print(permuted)

        PSIS_XY(polys=polys, m=m, permuted=permuted, permute_list=permute_list, levels=levels,
                recovered_path=recovered_path)

    # poly_names = [[],[],["F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9"],
    #               ["F10", "F11", "F12", "F13"],
    #               ["F14", "F15", "F16"]]
    # permuted_actions = [True]
    # for lvs in range(2, 3):
    #     lvs_poly_names = poly_names[lvs]
    #     for name in lvs_poly_names:
    #         if name != "F1":
    #             continue
    #         for permuted in permuted_actions:
    #             polys, m, permute_list = config_reader(name, lvs)
    #             recovered_path = RECOVERED_PREFIX + "/" + name + "/"
    #             if not permuted:
    #                 recovered_path = RECOVERED_PREFIX + "/" + name + "_w/"
    #             # print(permuted)
    #
    #             PSIS_XY(polys=polys, m=m, permuted=permuted, permute_list=permute_list, levels=lvs,
    #                     recovered_path=recovered_path)
    # invs = [[{}, {}, {}], [{}, {0: 1}, {0: 1}], [{}, {20: 1, 1: 1}, {23: 1, 21: 1, 20: 1, 4: 1, 1: 1, 0: 1}], [{}, {19: 1, 0: 1}, {23: 1, 22: 1, 21: 1, 19: 1, 4: 1, 3: 1, 1: 1}], [{}, {20: 1, 19: 1, 18: 1, 17: 1, 16: 1, 15: 1, 14: 1, 13: 1, 12: 1, 11: 1, 10: 1, 9: 1, 8: 1, 7: 1, 6: 1, 5: 1, 4: 1, 3: 1, 2: 1}, {23: 1, 22: 1, 20: 1, 19: 1, 18: 1, 17: 1, 16: 1, 15: 1, 14: 1, 13: 1, 12: 1, 11: 1, 10: 1, 9: 1, 8: 1, 7: 1, 6: 1, 5: 1, 1: 1}], [{}, {20: 1, 18: 1, 16: 1, 14: 1, 12: 1, 10: 1, 8: 1, 6: 1, 4: 1, 2: 1, 1: 1, 0: 1}, {23: 1, 20: 1, 18: 1, 16: 1, 14: 1, 12: 1, 10: 1, 8: 1, 6: 1, 1: 1, 0: 1}], [{}, {19: 1, 18: 1, 17: 1, 16: 1, 15: 1, 14: 1, 13: 1, 12: 1, 11: 1, 10: 1, 9: 1, 8: 1, 7: 1, 6: 1, 5: 1, 4: 1, 3: 1, 2: 1, 1: 1}, {22: 1, 21: 1, 19: 1, 18: 1, 17: 1, 16: 1, 15: 1, 14: 1, 13: 1, 12: 1, 11: 1, 10: 1, 9: 1, 8: 1, 7: 1, 6: 1, 5: 1, 4: 1, 0: 1}], [{}, {20: 1, 19: 1, 17: 1, 16: 1, 14: 1, 13: 1, 11: 1, 10: 1, 8: 1, 7: 1, 5: 1, 4: 1, 2: 1, 0: 1}, {22: 1, 21: 1, 20: 1, 19: 1, 17: 1, 16: 1, 14: 1, 13: 1, 11: 1, 10: 1, 8: 1, 7: 1, 5: 1, 4: 1, 3: 1, 1: 1}]]
    # levels = 2
    #
    # pixels_share1 = open_image("img/share/p1_share2.png")
    # pixels_share2 = open_image("img/share/p1_share3.png")
    #
    # pss = [pixels_share1, pixels_share2]
    #
    # mod_poly = {24: 1}
    # # 这里不包括完全恢复
    # rr = remove_matching_terms(get_pattern(3), mod_poly)
    # print("unmodified bits =", list(rr.keys()))
    # # rec = share_recover_2(pss, ids, mod_poly, invs, rr.keys(), lv,
    # #                       m=m, permuted=permuted, permute_list=permute_list)
    #
    # rec = share_recover_3(pss, [2,3], mod_poly, invs, rr.keys(), 2,
    #                       m=3, permuted=, permute_list=[0, 1, 2, 3, 4, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 5, 6, 16])
    # print(rec)
    # rec = pixels_1D_to_2D(rec, 512, 512)
    #
    # rec_np = np.array(rec)
    #
    # # 解密阶段
    # rec_uncrambled = None
    # if SCRAMBLED_OPTION == 'RANDOM':
    #     rec_uncrambled = unscramble_image_with_seed(rec_np, SEED)
    # if SCRAMBLED_OPTION == 'ARNOLD':
    #     rec_uncrambled = arnold_unscramble(rec_np, ARNOLD_ITERATIONS)
    #
    # res_path = 'img/out/recs.png'
    # save_share(res_path, rec_uncrambled)
