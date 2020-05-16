#! /usr/bin/env python
# -*- coding: utf-8 -*-
# coding=utf-8
import math
from math import ceil

count_of_worchop = 16

# transformers = [100, 250, 400, 630, 1000, 1600, 2500, 2*100, 2*250, 2*400, 2*630, 2*1000, 2*1600, 2*2500, 'РП', '2xРП']

# categories_arr = [2, 2, 3, 3, 2, 2, 2, 3, 2, 2, 3, 1, 1, 3, 3, 3 ,3 ,3 ]

# номинальная мощность
p_nom = [4750, 7600, 1460, 1390, 1540, 2200, 120,
         950, 8200, 1300, 1130, 380, 120, 300, 500, 300]

# коэффициент спроса ks
ks = [0.72, 0.89, 0.4, 0.69, 0.8, 0.5, 0.45, 0.5,
      0.5, 0.81, 0.75, 0.66, 0.84, 0.84, 0.91, 0.66]

price_of_electricity = 9
koef_amort = 0.03
working_hours = 2500
g1km = 1.48

# расчетная активная нагрузка
p_r = [0]*count_of_worchop
for i in range(count_of_worchop):
    p_r[i] = ceil(p_nom[i] * ks[i])
print(p_r, 'расчетная активная нагрузка')

# cosfi
cos_fi = [0.86, 0.6, 0.86, 0.86, 0.8, 0.8, 0.9,
          0.8, 0.86, 0.6, 0.8, 0.6, 0.6, 0.6, 0.86, 0.86]

# расчетная реактивная нагрузка
q_r = [0]*count_of_worchop
for i in range(count_of_worchop):
    q_r[i] = p_r[i] * math.tan(math.acos(cos_fi[i]))
# print(sum(q_r), 'расчетная реактивная нагрузка')

# расчетная полная мощность
s_r = [0]*count_of_worchop
for i in range(count_of_worchop):
    s_r[i] = math.sqrt(p_r[i]**2 + q_r[i]**2)
# print(sum(s_r), 's_r')
# print(s_r, 's_r')

# мощность компенсирующих устройств
tan_fz = 0.33
tan_fn = sum(q_r) / sum(p_r)
print(tan_fn, 'tan')
q_ku = sum(p_r) * (tan_fn - tan_fz)
# print(q_ku, 'q_ku')

# потери мощности в компенсирующем устройстве
d_p_ku = 0.002 * q_ku
# print(d_p_ku, 'p_ku')

# расчѐтная мощность на шинах ГПП
k_pm = 0.9
# print(sum(p_r), 'p_r')
# print(sum(q_r), 'q_r')
p_rgpp = (sum(p_r) + d_p_ku) * k_pm
q_rgpp = (sum(q_r) - q_ku) * k_pm

# полная расчѐтная мощность на шинах ГПП
s_rgpp = math.sqrt(p_rgpp**2 + q_rgpp**2)
# print(s_rgpp, 's_rgpp')

# потери мощности в трансформаторе ГПП
d_r_tgpp = 0.02 * s_rgpp
# print(d_r_tgpp, 'd_r_tgpp')
d_q_tpgg = 0.01 * s_rgpp

# полная расчетная нагрузка с учѐтом потерь мощности
s__rgpp = math.sqrt((p_rgpp + d_r_tgpp)**2 + (q_rgpp + d_q_tpgg)**2)
# print(s__rgpp, 's__rgpp')

# Выбор высоковольтных потребителей

# Выбор напряжения питающей и распределительной сети. Выбор напряжения питающей линии

l = 0.6  # длина линии

u1 = 3 * math.sqrt(sum(s_r) * 10**-3) + 0.5 * l
u2 = 4.34 * math.sqrt(l + 16 * sum(p_r) * 10**-3)
u3 = 16 * math.sqrt(math.sqrt(sum(p_r) * 10**-3 * l))
u4 = 17 * math.sqrt((l/16) + sum(p_r) * 10**-3)

u = (u1 + u2 + u3 + u4) / 4
# print(u, 'Напряжение срееднее арифм')

# Определение типа, количества и мощности цеховых трансформаторных подстанций

# print(p_r[0], 'perviy zeh')
# print(q_r[0], 'perviy zeh')

tan_fi_z = 0.33

# Определим мощность, тип, и количество компенсирующих устройств

tan_fi1 = [0]*count_of_worchop
for i in range(count_of_worchop):
    tan_fi1[i] = q_r[i] / p_r[i]

q_ku_1 = [0]*count_of_worchop
# q_ku_choosen  = [0]*count_of_worchop
for i in range(count_of_worchop):
    q_ku_1[i] = ceil(p_r[i]*(tan_fi1[i] - tan_fi_z))

# print(q_ku_1, 'q_ku_1')

# комплексную реактивную мощность после компенсации
transformers_table = [2500, 2500, 1000, 1600, 1000,
                      1000, 1, 630, 1600, 1600, 630, 2, 630, 400, 1000, 1]
transformers_count = [2, 4, 1, 1, 2, 2, 0, 1, 4, 1, 2, 0, 1, 1, 1, 0]
# print(transformers_count)
# print(len(transformers_count), 'transformers count')
q_ku_chosen = [0]*count_of_worchop
for i in range(count_of_worchop):
    if q_ku_1[i] <= 25 and transformers_count[i] != 0:
        q_ku_chosen[i] = 25
    elif transformers_count[i] != 0:
        q_ku_chosen[i] = math.ceil(
            ((q_ku_1[i]) / (transformers_count[i]*25)))*25
# for i in range(count_of_worchop):
#     for j in transformers_count:
#             if q_ku_1[i]<=25:
#                 q_ku_chosen[i]=25
#             else:
#                 q_ku_chosen[i] =ceil (((q_ku_1[i]) // (transformers_count[j]*25)))*25

# print(q_ku_chosen, 'q_ku_chosen')

# комплексная реактивная мощность после компенсации
q__p = [0]*count_of_worchop
for i in range(count_of_worchop):
    q__p[i] = ceil(q_r[i] - (transformers_count[i] * q_ku_chosen[i]))

# print(q__p, 'комплексная реактивная мощность после компенсации')

# расчетная активная мощность
p__r = [0]*count_of_worchop
for i in range(count_of_worchop):
    if i == 3:
        p__r[3] = ceil(p_r[3] + p_r[6])
    elif i == 12:
        p__r[12] = ceil(p_r[12] + p_r[11])
    elif i == 14:
        p__r[14] = ceil(p_r[14] + p_r[15])
    else:
        p__r[i] = ceil(p_r[i])
print(p__r, 'p__r')

# расчетная реактивная мощность
q__r = [0]*count_of_worchop
for i in range(count_of_worchop):
    if i == 3:
        q__r[3] = ceil(q_r[3] + q_r[6])
    elif i == 12:
        q__r[12] = ceil(q_r[12] + q_r[11])
    elif i == 14:
        q__r[14] = ceil(q_r[14] + q_r[15])
    else:
        q__r[i] = ceil(q_r[i])
print(q__r, 'q__r')
# Полная расчѐтная нагрузка трансформаторов
s_p = [0]*count_of_worchop
for i in range(count_of_worchop):
    s_p[i] = ceil(math.sqrt(p__r[i]**2 + q__p[i]**2))
# print(s_p, 'Полная расчѐтная нагрузка трансформаторов')

# Выбор кабельных линий распределительной сети
p_hh = [
    2.8,  # 2500
    4.8,  # тсз 2500
    1.5,  # 1000
    2.6,  # 1600
    1.5,  # 1000
    1.5,  # 1000
    0,
    1.05,  # 630
    2.6,  # 1600
    2.6,  # 1600
    1.05,  # 630
    0,
    1.05,  # 630
    0.83,  # 400
    1.5,  # 1000
    0
]

p_kz = [
    24,  # 2500
    12,  # тсз 2500
    10.5,  # 1000
    17,  # 1600
    10.5,  # 1000
    10.5,  # 1000
    0,
    6.5,  # 630
    17,  # 1600
    17,  # 1600
    6.5,  # 630
    0,
    6.5,  # 630
    4.6,  # 400
    10.5,  # 1000
    0
]

u_kz = [
    6,  # 2500
    6,  # тсз 2500
    5.5,  # 1000
    6,  # 1600
    5.5,  # 1000
    5.5,  # 1000
    0,
    5.5,  # 630
    6,  # 1600
    6,  # 1600
    5.5,  # 630
    0,
    5.5,  # 630
    4.5,  # 400
    5.5,  # 1000
    0
]

i_xx = [
    0.8,  # 2500
    0.75,  # тсз 2500
    1.2,  # 1000
    1,  # 1600
    1.2,  # 1000
    1.2,  # 1000
    0,
    1.8,  # 630
    1,  # 1600
    1,  # 1600
    1.8,  # 630
    0,
    1.8,  # 630
    2,  # 400
    1.2,  # 1000
    0
]

kzn = [0] * count_of_worchop
for i in range(count_of_worchop):
    if transformers_table[i] >= 10:
        kzn[i] = s_p[i]/(transformers_count[i]*transformers_table[i])
    else:
        pass
    # kzn[i] = s_p[i]/(transformers_count[i]*transformers_table[i])
print(kzn, "kzn")

# Выбор кабельных линий распределительной сети

# Определим потери активной мощности ТП
delta_p_tpi = [0]*count_of_worchop
for i in range(count_of_worchop):
    delta_p_tpi[i] = p_hh[i] + (p_kz[i] * kzn[i]**2)
print(delta_p_tpi, "delta_p_tpi")

# Определим потери реактивной мощности ТП

q_hh = [0] * count_of_worchop
for i in range(count_of_worchop):
    if transformers_table[i] >= 10:
        q_hh[i] = (i_xx[i]/100) * transformers_table[i]
# print(q_hh, "q_hh")

q_kz = [0] * count_of_worchop
for i in range(count_of_worchop):
    if transformers_table[i] >= 10:
        q_kz[i] = (u_kz[i]/100) * transformers_table[i]

delta_q_tpi = [0]*count_of_worchop
for i in range(count_of_worchop):
    delta_q_tpi[i] = q_hh[i] + (q_kz[i] * kzn[i]**2)

print(delta_q_tpi, "delta_q_tpi")
# Определяем расчѐтную нагрузку линии ГПП - ТП

s_r_a_gpp_tpi = [0] * count_of_worchop
for i in range(count_of_worchop):
    s_r_a_gpp_tpi[i] = round(
        math.sqrt((p__r[i] + delta_p_tpi[i])**2 + (q__r[i] + delta_q_tpi[i])**2), 3)
print(s_r_a_gpp_tpi, 's_r_a_gpp_tpi')

s_r_n_gpp_tpi = [0] * count_of_worchop
for i in range(count_of_worchop):
    if transformers_count[i] != 0:
        s_r_n_gpp_tpi[i] = round(math.sqrt(
            ((p__r[i]/transformers_count[i]) + delta_p_tpi[i])**2 + ((q__r[i]/transformers_count[i]) + delta_q_tpi[i])**2), 3)
    else:
        pass
print(s_r_n_gpp_tpi, 's_r_n_gpp_tpi')

# Определяем расчѐтную токовую нагрузку линии ГПП - ТП

i_r_n_gpp_tpi = [0] * count_of_worchop
for i in range(count_of_worchop):
    i_r_n_gpp_tpi[i] = round(s_r_n_gpp_tpi[i]/(math.sqrt(3)*10), 3)
print(i_r_n_gpp_tpi, 'i_r_n_gpp_tpi')
i_r_a_gpp_tpi = [0] * count_of_worchop
for i in range(count_of_worchop):
    i_r_a_gpp_tpi[i] = round(s_r_a_gpp_tpi[i]/(math.sqrt(3)*10), 3)
print(i_r_a_gpp_tpi, 'i_r_a_gpp_tpi')

# Определяем допустимую токовую нагрузку линии ГПП - ТП

# Если количество кабелей 1, то коэффициент 1.. Если 2, то 0.9. 4 - 0,8
k1 = [0.9, 0.8, 1, 1, 0.9, 0.9, 1, 1, 0.8, 1, 0.9, 1, 1, 1, 1, 1]
k2 = 1.15
i_dop_n_gpp_tpi = [0] * count_of_worchop
for i in range(count_of_worchop):
    i_dop_n_gpp_tpi[i] = i_r_n_gpp_tpi[i]/k1[i]

i_dop_a_gpp_tpi = [0] * count_of_worchop
for i in range(count_of_worchop):
    i_dop_a_gpp_tpi[i] = round(i_r_a_gpp_tpi[i]/(k1[i]*k2), 3)
print(i_dop_a_gpp_tpi, 'i_dop_a_gpp_tpi')
# Выбираем сечение линии по нагреву

# Длительно-допустимые токовые нагрузки кабелей
i_dop_kab = [
    74,  # 16мм2
    91,  # 25 мм2
    110,  # 35
    134,  # 50
    162,  # 70
    192,  # 95
    218,  # 120
    246,  # 150
    275,  # 185 мм2
    314,  # 240
]


class Cable:

    # i_dop
    # i_dop_kab_price
    # sechenie
    # l_delta_u1percent
    def __init__(self,
                 i_dop,
                 i_dop_kab_price,
                 sechenie,
                 l_delta_u1percent):
        self.i_dop = i_dop
        self.i_dop_kab_price = i_dop_kab_price
        self.sechenie = sechenie
        self.l_delta_u1percent


# class Asb(Cable):
#     def __init__(self,
#                  i_dop,
#                  i_dop_kab_price,
#                  sechenie,
#                  l_delta_u1percent,
#                  s_asb):
#         super().__init__(self,
#                          i_dop,
#                          i_dop_kab_price,
#                          sechenie,
#                          l_delta_u1percent,
#                          ):
#                     self.s_asb = s_asb


class Aasv(Cable):
    pass


asb_arr = []


# def init_asb(parameter_list):
#     for i in range(len(i_dop_kab_aasb)):


i_dop_kab_price = [
    165,
    175,
    203,
    242,
    299,
    366,
    381,
    500,
    588,
    671,
]

i_dop_kab_aasb = [
    79,  # 16
    102,  # 25,
    126,  # 35,
    153,  # 50,
    184,  # 70,
    219,  # 95,
    248,  # 120,
    281,  # 150,
    314,  # 185,
    359  # 240
]

i_dop_kab_aasb_price = [
    385,
    331,
    371,
    435,
    504,
    601,
    676,
    546,
    1024,
    1091
]

s_asb = [
    51,  # кВт
    67,
    82,
    100,
    120,
    143,
    163,
    184,
    206,
    236
]

l_delta_u1percent = [
    95.1,  # 16
    148.6,  # 25,
    208.1,  # 35,
    297.3,  # 50,
    416.2,  # 70,
    564.9,  # 95,
    713.5,  # 120,
    892,  # 150,
    1100,  # 185,
    1427  # 240
]

длины_линий = [
    998,
    934,
    770,
    750,
    434,
    223,
    619,
    474,
    231,
    620,
    447,
    201,
    447,
    497,
    105,
    179
]

l_delta_u1percent_chosen = [0] * count_of_worchop

chosen_cables = [0] * count_of_worchop

for i in range(count_of_worchop):
    if i_dop_a_gpp_tpi[i] <= 74 and transformers_table[i] > 0:
        chosen_cables[i] = 16
        if chosen_cables[i] == 16:
            l_delta_u1percent_chosen[i] = l_delta_u1percent[0]

    elif i_dop_a_gpp_tpi[i] > 74 and i_dop_a_gpp_tpi[i] <= 91 and transformers_table[i] > 0:
        chosen_cables[i] = 25
        if chosen_cables[i] == 25:
            l_delta_u1percent_chosen[i] = l_delta_u1percent[1]
    elif i_dop_a_gpp_tpi[i] > 91 and i_dop_a_gpp_tpi[i] <= 110 and transformers_table[i] > 0:
        chosen_cables[i] = 35
        if chosen_cables[i] == 35:
            l_delta_u1percent_chosen[i] = l_delta_u1percent[2]
    elif i_dop_a_gpp_tpi[i] > 110 and i_dop_a_gpp_tpi[i] <= 134 and transformers_table[i] > 0:
        chosen_cables[i] = 50
        if chosen_cables[i] == 50:
            l_delta_u1percent_chosen[i] = l_delta_u1percent[3]
    elif i_dop_a_gpp_tpi[i] > 134 and i_dop_a_gpp_tpi[i] <= 162 and transformers_table[i] > 0:
        chosen_cables[i] = 70
        if chosen_cables[i] == 70:
            l_delta_u1percent_chosen[i] = l_delta_u1percent[4]
    elif i_dop_a_gpp_tpi[i] > 162 and i_dop_a_gpp_tpi[i] <= 192 and transformers_table[i] > 0:
        chosen_cables[i] = 95
        if chosen_cables[i] == 95:
            l_delta_u1percent_chosen[i] = l_delta_u1percent[5]
    elif i_dop_a_gpp_tpi[i] > 192 and i_dop_a_gpp_tpi[i] <= 218 and transformers_table[i] > 0:
        chosen_cables[i] = 120
        if chosen_cables[i] == 120:
            l_delta_u1percent_chosen[i] = l_delta_u1percent[6]
    elif i_dop_a_gpp_tpi[i] > 218 and i_dop_a_gpp_tpi[i] <= 246 and transformers_table[i] > 0:
        chosen_cables[i] = 150
        if chosen_cables[i] == 150:
            l_delta_u1percent_chosen[i] = l_delta_u1percent[7]
    elif i_dop_a_gpp_tpi[i] > 246 and i_dop_a_gpp_tpi[i] <= 275 and transformers_table[i] > 0:
        chosen_cables[i] = 185
        if chosen_cables[i] == 185:
            l_delta_u1percent_chosen[i] = l_delta_u1percent[8]
    elif i_dop_a_gpp_tpi[i] > 275:
        chosen_cables[i] = 240
        if chosen_cables[i] == 240:
            l_delta_u1percent_chosen[i] = l_delta_u1percent[9]

    elif i_dop_a_gpp_tpi[i] <= 79 and s_r_a_gpp_tpi[i] <= 51 and transformers_table[i] == 0:
        chosen_cables[i] = 16
    elif i_dop_a_gpp_tpi[i] > 79 and i_dop_a_gpp_tpi[i] <= 102 or (s_r_a_gpp_tpi[i] > 51 and s_r_a_gpp_tpi[i] <= 67) and transformers_table[i] == 0:
        chosen_cables[i] = 25
    elif i_dop_a_gpp_tpi[i] > 102 and i_dop_a_gpp_tpi[i] <= 126 or (s_r_a_gpp_tpi[i] > 67 and s_r_a_gpp_tpi[i] <= 82) and transformers_table[i] == 0:
        chosen_cables[i] = 35
    elif i_dop_a_gpp_tpi[i] > 126 and i_dop_a_gpp_tpi[i] <= 153 or (s_r_a_gpp_tpi[i] > 82 and s_r_a_gpp_tpi[i] <= 100) and transformers_table[i] == 0:
        chosen_cables[i] = 50
    elif i_dop_a_gpp_tpi[i] > 153 and i_dop_a_gpp_tpi[i] <= 184 or (s_r_a_gpp_tpi[i] > 100 and s_r_a_gpp_tpi[i] <= 120) and transformers_table[i] == 0:
        chosen_cables[i] = 70
    elif i_dop_a_gpp_tpi[i] > 184 and i_dop_a_gpp_tpi[i] <= 219 or (s_r_a_gpp_tpi[i] > 120 and s_r_a_gpp_tpi[i] <= 143) and transformers_table[i] == 0:
        chosen_cables[i] = 95
    elif i_dop_a_gpp_tpi[i] > 219 and i_dop_a_gpp_tpi[i] <= 248 or (s_r_a_gpp_tpi[i] > 143 and s_r_a_gpp_tpi[i] <= 163) and transformers_table[i] == 0:
        chosen_cables[i] = 120
    elif i_dop_a_gpp_tpi[i] > 248 and i_dop_a_gpp_tpi[i] <= 281 or (s_r_a_gpp_tpi[i] > 163 and s_r_a_gpp_tpi[i] <= 184) and transformers_table[i] == 0:
        chosen_cables[i] = 150
    elif i_dop_a_gpp_tpi[i] > 281 and i_dop_a_gpp_tpi[i] <= 314 or (s_r_a_gpp_tpi[i] > 184 and s_r_a_gpp_tpi[i] <= 206) and transformers_table[i] == 0:
        chosen_cables[i] = 185
    elif i_dop_a_gpp_tpi[i] > 314 or s_r_a_gpp_tpi[i] > 206 and transformers_table[i] == 0:
        chosen_cables[i] = 240

i_dop_of_chosen_cables = [0] * count_of_worchop
for i in range(count_of_worchop):
    if chosen_cables[i] == 16 and transformers_count[i] > 0:
        i_dop_of_chosen_cables[i] = 74
    elif chosen_cables[i] == 25 and transformers_count[i] > 0:
        i_dop_of_chosen_cables[i] = 91
    elif chosen_cables[i] == 35 and transformers_count[i] > 0:
        i_dop_of_chosen_cables[i] = 110
    elif chosen_cables[i] == 50 and transformers_count[i] > 0:
        i_dop_of_chosen_cables[i] = 134
    elif chosen_cables[i] == 70 and transformers_count[i] > 0:
        i_dop_of_chosen_cables[i] = 162
    elif chosen_cables[i] == 95 and transformers_count[i] > 0:
        i_dop_of_chosen_cables[i] = 192
    elif chosen_cables[i] == 120 and transformers_count[i] > 0:
        i_dop_of_chosen_cables[i] = 218
    elif chosen_cables[i] == 150 and transformers_count[i] > 0:
        i_dop_of_chosen_cables[i] = 246
    elif chosen_cables[i] == 185 and transformers_count[i] > 0:
        i_dop_of_chosen_cables[i] = 275
    elif chosen_cables[i] == 240 and transformers_count[i] > 0:
        i_dop_of_chosen_cables[i] = 314

    elif chosen_cables[i] == 16 and transformers_count[i] == 0:
        i_dop_of_chosen_cables[i] = 79
    elif chosen_cables[i] == 25 and transformers_count[i] == 0:
        i_dop_of_chosen_cables[i] = 102
    elif chosen_cables[i] == 35 and transformers_count[i] == 0:
        i_dop_of_chosen_cables[i] = 126
    elif chosen_cables[i] == 50 and transformers_count[i] == 0:
        i_dop_of_chosen_cables[i] = 153
    elif chosen_cables[i] == 70 and transformers_count[i] == 0:
        i_dop_of_chosen_cables[i] = 184
    elif chosen_cables[i] == 95 and transformers_count[i] == 0:
        i_dop_of_chosen_cables[i] = 219
    elif chosen_cables[i] == 120 and transformers_count[i] == 0:
        i_dop_of_chosen_cables[i] = 248
    elif chosen_cables[i] == 150 and transformers_count[i] == 0:
        i_dop_of_chosen_cables[i] = 281
    elif chosen_cables[i] == 185 and transformers_count[i] == 0:
        i_dop_of_chosen_cables[i] = 314
    elif chosen_cables[i] == 240 and transformers_count[i] == 0:
        i_dop_of_chosen_cables[i] = 359


print(chosen_cables, 'chosen_cables')
print(l_delta_u1percent_chosen, 'l_delta_u1percent_chosen')
l_dop = [0] * count_of_worchop
for i in range(count_of_worchop):
    l_dop[i] = round(l_delta_u1percent_chosen[i] *
                     10 * (chosen_cables[i]/i_dop_a_gpp_tpi[i]), 3)

s_e = [0] * count_of_worchop
for i in range(count_of_worchop):
    s_e[i] = round(i_r_n_gpp_tpi[i]/1.4, 3)

print(s_e, 's_e')
print(l_dop, 'l_dop')
print(s_asb, 's_asb')

prices_of_chosen_cables = [0] * count_of_worchop
for i in range(count_of_worchop):
    if chosen_cables[i] == 16 and transformers_count[i] > 0:
        prices_of_chosen_cables[i] = i_dop_kab_price[0]
    elif chosen_cables[i] == 25 and transformers_count[i] > 0:
        prices_of_chosen_cables[i] = i_dop_kab_price[1]
    elif chosen_cables[i] == 35 and transformers_count[i] > 0:
        prices_of_chosen_cables[i] = i_dop_kab_price[2]
    elif chosen_cables[i] == 50 and transformers_count[i] > 0:
        prices_of_chosen_cables[i] = i_dop_kab_price[3]
    elif chosen_cables[i] == 70 and transformers_count[i] > 0:
        prices_of_chosen_cables[i] = i_dop_kab_price[4]
    elif chosen_cables[i] == 95 and transformers_count[i] > 0:
        prices_of_chosen_cables[i] = i_dop_kab_price[5]
    elif chosen_cables[i] == 120 and transformers_count[i] > 0:
        prices_of_chosen_cables[i] = i_dop_kab_price[6]
    elif chosen_cables[i] == 150 and transformers_count[i] > 0:
        prices_of_chosen_cables[i] = i_dop_kab_price[7]
    elif chosen_cables[i] == 185 and transformers_count[i] > 0:
        prices_of_chosen_cables[i] = i_dop_kab_price[8]
    elif chosen_cables[i] == 240 and transformers_count[i] > 0:
        prices_of_chosen_cables[i] = i_dop_kab_price[9]

    elif chosen_cables[i] == 16 and transformers_count[i] == 0:
        prices_of_chosen_cables[i] = i_dop_kab_aasb_price[0]
    elif chosen_cables[i] == 25 and transformers_count[i] == 0:
        prices_of_chosen_cables[i] = i_dop_kab_aasb_price[1]
    elif chosen_cables[i] == 35 and transformers_count[i] == 0:
        prices_of_chosen_cables[i] = i_dop_kab_aasb_price[2]
    elif chosen_cables[i] == 50 and transformers_count[i] == 0:
        prices_of_chosen_cables[i] = i_dop_kab_aasb_price[3]
    elif chosen_cables[i] == 70 and transformers_count[i] == 0:
        prices_of_chosen_cables[i] = i_dop_kab_aasb_price[4]
    elif chosen_cables[i] == 95 and transformers_count[i] == 0:
        prices_of_chosen_cables[i] = i_dop_kab_aasb_price[5]
    elif chosen_cables[i] == 120 and transformers_count[i] == 0:
        prices_of_chosen_cables[i] = i_dop_kab_aasb_price[6]
    elif chosen_cables[i] == 150 and transformers_count[i] == 0:
        prices_of_chosen_cables[i] = i_dop_kab_aasb_price[7]
    elif chosen_cables[i] == 185 and transformers_count[i] == 0:
        prices_of_chosen_cables[i] = i_dop_kab_aasb_price[8]
    elif chosen_cables[i] == 240 and transformers_count[i] == 0:
        prices_of_chosen_cables[i] = i_dop_kab_aasb_price[9]

print(prices_of_chosen_cables, 'prices_of_chosen_cables')


# Технико-экономический расчѐт кабельных линий

# Капитальные вложения на линию:
nk = 2  # колво
k_gpp_tpi = [0] * count_of_worchop
for i in range(count_of_worchop):
    k_gpp_tpi[i] = nk * длины_линий[i] * \
        pow(10, -3) * prices_of_chosen_cables[i]
print(k_gpp_tpi, 'k_gpp_tpi')

# расход цветного металла

g_gpp_tpi = [0] * count_of_worchop
for i in range(count_of_worchop):
    g_gpp_tpi[i] = nk * длины_линий[i] * pow(10, -3) * g1km

# Коэффициент загрузки и расчѐтная полная нагрузка рассматриваемой линии

kz = [0] * count_of_worchop
for i in range(count_of_worchop):
    kz[i] = i_r_n_gpp_tpi[i] / i_dop_of_chosen_cables[i]
print(kz, 'kz')

# Потери мощности рассматриваемой линии:
delta_p_gpp_tpi = [0] * count_of_worchop
for i in range(count_of_worchop):
    delta_p_gpp_tpi[i] = nk * длины_линий[i] * pow(10, -3) * kz[i]**2

# Потери электроэнергии:
delta_e_gpp_tpi = [0] * count_of_worchop
for i in range(count_of_worchop):
    delta_e_gpp_tpi[i] = round(delta_p_gpp_tpi[i] * working_hours, 3)
print(delta_e_gpp_tpi, 'delta_e_gpp_tpi')

# Стоимость потерь электроэнергии:
sp_gpp_tpi = [0] * count_of_worchop
for i in range(count_of_worchop):
    sp_gpp_tpi[i] = delta_e_gpp_tpi[i] * price_of_electricity * pow(10, -3)
print(sp_gpp_tpi, 'sp_gpp_tpi')

# Стоимость амортизационных отчислений на рассматриваемую линию:
sa_gpp_tpi = [0] * count_of_worchop
for i in range(count_of_worchop):
    sa_gpp_tpi[i] = round(koef_amort * k_gpp_tpi[i], 2)
print(sa_gpp_tpi, 'sa_gpp_tpi')
