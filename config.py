from config_poly_group import *

RECOVER_LEVEL = 1

# 配置需要的多项式组合
LEVELS = 2  # 要实现几级恢复
config_name = "F1"  # 使用的多项式组
config_group = ALL_GROUP.get(LEVELS)

cfg = config_group.get(config_name)
CONFIG_G1 = cfg.get("G1")
CONFIG_G2 = cfg.get("G2")
CONFIG_G3 = cfg.get("G3")
CONFIG_G4 = cfg.get("G4")
M = cfg.get("M")  # 一个参数内处理几个像素

K = 2  # K和每次要提供的id个数是一致的 TODO 目前只支持K=2和K=3

# 预先置乱，默认OPEN
PERMUTED = True
PERMUTE_LIST = cfg.get("permute_list")

# 对modified bits的设置，置0/置1/置随机0、1
MODIFIED_PIXEL_SET = {
    0: 'BLACK',
    1: 'WHITE',
    2: 'RANDOM',
    3: 'RANDOM_OVERALL'
}
MODIFIED_PIXEL_OPTION = 0  # 默认0


SCRAMBLED_TYPES = {
    0: 'RANDOM',
    1: 'ARNOLD'
}

ARNOLD_ITERATIONS = 5
SEED = 42

SHARE_SUFFIX = '.png'
RECOVERD_SUFFIX = '.png'
SCRAMBLED_OPTION = SCRAMBLED_TYPES.get(0)  # 用随机置乱

MODE_SMALL_TEST = True


FILE_FOLDER_NAME = config_name
if not PERMUTED:
    FILE_FOLDER_NAME = FILE_FOLDER_NAME + "_w"

DATASETS = {
    0: "dataset_480_p500",
    1: "dataset_512_p49",
    2: "dataset_512_p100"
}

LAB_ROOT_PATH = "E:/YURI_research"
DATASET = DATASETS.get(2)

if MODE_SMALL_TEST:
    TEST_IMG_PATH = 'img/test_img/boss/'
    SHARE_PATH = 'img/share/'
    SUB_SHARE_PATH = 'img/share/'
    RECOVERED_PATH = 'img/out/'
    RECOVERED_PREFIX = 'img/out'
    ANALYSED_PATH = 'img/out/analysed/'
    REF_PATH = 'reference/'
else:
    TEST_IMG_PATH = f'{LAB_ROOT_PATH}/experiment/secret/{DATASET}'
    SHARE_PATH = f'{LAB_ROOT_PATH}/experiment/result/share/'
    SUB_SHARE_PATH = f'{LAB_ROOT_PATH}/experiment/result/share/'
    RECOVERED_PATH = f'{LAB_ROOT_PATH}/experiment/result/{DATASET}/{FILE_FOLDER_NAME}/'
    RECOVERED_PREFIX = f'{LAB_ROOT_PATH}/experiment/result/{DATASET}'
    ANALYSED_PATH = f'{LAB_ROOT_PATH}/experiment/result/analysed/'
    REF_PATH = f'E:/YURI_research/experiment/reference/{DATASET}/'
