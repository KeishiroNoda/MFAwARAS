import struct
import random
import numpy
import math

# 64ビット（8バイト）浮動小数点数
#    実数値を (-1)^{符号} x (仮数) x 2^{指数}の形式で表現．ただし仮数は1以上2未満に取る
# 63 (MSB): 符号．非負なら0, 負なら1
# 62-52: 指数表現部．（指数 + 1023)を11ビットで表現する
#    1.0 = 1.0 x 2^{0} なら 0+1023 = 1023 = (011 1111 1111)_b = (3 ff)_x
#    0.5 = 1.0 x 2^{-1} なら -1+1023 = 1022 = (011 1111 1110)_b = (3 fe)_x
#    2.0 = 1.0 x 2^{1} なら 1+1023 = 1024 = (100 0000 0000)_b = (4 00)_x
#    (0 00)_xは0, (7ff)_xは無限大を示す特殊値
# 51-0 (LSB)：仮数表現部．仮数（1以上2未満なので，必ず 1.xxxx...の形式）の小数部分を52ビットで表現する
#    1.5 = (1.1)_b -> (1000 ...0)_b = (80 00...00)_x
#    1.25 = (1.01)_b -> (0100 ...0)_b = (40 00...00)_x
#    1.125 = (1.001)_b -> (0010 ...0)_b = (20 00...00)_x
#    1.75 = (1.11)_b -> (1100 ...0)_b = (c0 00...00)_x

# 浮動小数点数を，符号（1ビット），指数表現部（11ビット），仮数表現部（52ビット）の3つの部分に分解する
# 戻り値は整数のタプル（順に符号，指数表現部，仮数表現部の整数表現）
def double_partition(f):
    # 浮動小数点を，符号なしlong long intに変換
    expanded = struct.unpack('>Q', struct.pack('>d', f))[0]
    frac_part = expanded & ((1 << 52) - 1) # 下位52ビットを取得し frac_partに
    expanded = expanded >> 52 # 上位12ビットを下位に降ろす
    exp_part = expanded & ((1 << 11) - 1) # 降ろした12ビットの下位11ビットを取得し exp_partに
    expanded = expanded >> 11 # 最上位ビットを下位に降ろす
    sig_part = expanded & 0b1 # 降ろした最上位ビットを sig_partに
    return sig_part, exp_part, frac_part

# double_partitionの結果をhexまたはbin表示
def double_partition_hex(f):
    sig_part, exp_part, frac_part = double_partition(f)
    return format(sig_part, '01b'), format(exp_part, '03x'), format(frac_part, '011x')

def double_partition_bin(f):
    sig_part, exp_part, frac_part = double_partition(f)
    return format(sig_part, '01b'), format(exp_part, '011b'), format(frac_part, '052b')

# 浮動小数点数（0より大きく 1未満であることを仮定）の小数部を，52ビット長のビット列（文字列）として返す
def double_to_binstr(f):
    sig_part, exp_part, frac_part = double_partition(f)
    frac_str = ''
    for i in range(1022 - exp_part):
        frac_str = frac_str + '0'
    frac_str = frac_str + '1' + format(frac_part, '052b')
    return frac_str[0:52]

# 52ビット長のビット列を小数部に持つ浮動小数点数を返す
def binstr_to_double(frac_str):
    frac_part = int(frac_str, 2)
    f = frac_part / (2 ** 52)
    return f

# 真値に対しノイズを足し（あるいは引き），真値に「遠すぎず近すぎない」近似値を作成したい．
# ノイズの絶対値は，小数点以下52ビットまで表示する2進表現で
#   0.00...010...0（下位k-1ビットすべて0, kビット目に1）以上
#   0.00...011...1（下位k-1ビットすべて1, kビット目に1）以下とする
# この範囲内のノイズであれば，0以上1以下の真値に作用させても「丸められない」（ただし，真値が0, 1に非常に近いときは要注意）
# ノイズの下位k-1ビットはランダムに選ぶことができるため，ノイズパターンは2^(k-1)通り存在
# プラス・マイナスの符号をランダムに与えると，全部で2^k通りのバリエーションを作り出すことができる
# このノイズを「k-ビット強度を持つノイズ」と呼ぶ

# 以下は，絶対値が最大となるノイズ，最小となるノイズ，k-ビット強度を持つランダムノイズを生成する関数
# いずれも，最初に整数としてノイズパターンを生成し，その整数を割ることで，小数点以下の部分にノイズパターンを持っていく
def gen_noise_absmax(k):
    base = 1 << k - 1
    noise = base + (1 << (k - 1)) - 1
    return (noise / base) / (2**(53 - (k - 1)))

def gen_noise_absmin(k):
    base = 1 << k - 1
    noise = base
    return (noise / base) / (2**(53 - (k - 1)))

def gen_noise(k):
    base = 1 << k - 1
    noise = base + random.randint(0, (1 << (k - 1)) - 1)
    if random.randint(0,1) == 1:
        noise = - noise
    return (noise / base) / (2**(53 - (k - 1)))

# 浮動小数点数からlenビットの鍵（小数部の上位lenビット）を取り出す
def get_key(val, len):
    return double_to_binstr(val)[0:len]

# 2つの文字列に対し，一致するプレフィックスの長さを返す
def count_match(p1, p2):
    length = min(len(p1), len(p2))
    i = 0
    while (i < length) and (p1[i] == p2[i]):
        i = i + 1
    return i

# ロジスティック写像の計算
def logistic_map(f):
    return 4.0 * f * (1 - f)

# ２つの初期値に対しロジスティック写像を繰り返し適用し，上位lenビットが一致しなくなる繰り返し回数を求める．
# 「繰り返し数が n」=「ロジスティック写像を n回適用すると不一致が発生」=「世代0から世代n-1までは上位lenビットが一致」
def match_until(a_init, b_init, len):
    a = a_init
    b = b_init
    count = 0
    # ロジスティック写像の計算結果が 0になると，その先ずっと 0から値が動かない．この場合，強制的にループから脱出する
    while (get_key(a, len) == get_key(b, len)) and (a != 0) and (b != 0):
        count = count + 1
        a = logistic_map(a)
        b = logistic_map(b)
    return count

# +absmax. +absmin, -absmin, -absmaxの４パターンのノイズをbaseに加え，何世代の同期が提供されるかを計算
# 同期世代数の最大と最小を返す
def ttl(base, noise_absmax, noise_absmin, bits_to_use):
    match_1 = match_until(base, base + noise_absmax, bits_to_use)
    match_2 = match_until(base, base + noise_absmin, bits_to_use)
    match_3 = match_until(base, base - noise_absmax, bits_to_use)
    match_4 = match_until(base, base - noise_absmin, bits_to_use)
    return max(match_1, match_2, match_3, match_4), min(match_1, match_2, match_3, match_4)

# メインの処理
def main():
    # ノイズのビット強度を指定
    strength = 16

    # 浮動小数点数の小数部の2進表示のうち，鍵として使用する上位ビット数
    # ここでは，strengthビットに相当するビット数としている
    bits_to_use = strength

    # 絶対値が最大・最小となるノイズを生成
    noise_absmax = gen_noise_absmax(strength)
    noise_absmin = gen_noise_absmin(strength)

    # 真値をランダムに選び，ノイズの影響（同期世代数の最大と最小）を評価する
    samples = 100
    for i in range(samples):
        base = random.uniform(0,1)
        ttl_max, ttl_min = ttl(base, noise_absmax, noise_absmin, bits_to_use)
        print(base, ttl_max, ttl_min)

if __name__ == "__main__":
    main()
