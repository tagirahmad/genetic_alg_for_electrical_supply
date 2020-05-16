# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-

from math import sqrt, ceil, acos, tan, floor
import itertools
from constants.constants import *
from constants.aasv_cable import *
from constants.asb_cable import *
from constants.circuit_breaker import *
from constants.transformer import *


class Цех:
    calculated_active_power = 0
    calculated_reactive_power = 0
    calculated_active_power_summ = 0
    calculated_reactive_power_сумма = 0
    calculated_full_power = 0
    tan_fi_n = 0
    power_of_compensating_devices = 0
    chosen_compens_device = 0
    complex_reactive_power_after_compensation = 0
    full_calculated_transformers_load = 0
    calculated_transformers_power = 0
    preselected_transformer = None
    selected_transformer = None
    transformer_load_coeff = 0
    active_power_loss_ts = None
    reactive_power_loss_ts = None
    permissible_voltage_loss = None
    permissible_economic_current_density_loss = None
    calculated_load_a_gpp_tp = None
    calculated_load_n_gpp_tp = None
    calculated_current_load_a_gpp_tp = None
    calculated_current_load_n_gpp_tp = None
    permissible_current_load_a_gpp_tp = None
    permissible_current_load_n_gpp_tp = None
    k1 = None
    k2 = 1.15
    cable_current = 0
    chosen_cable = None
    line_investments = 0
    color_material_consumption = 0
    load_coeff = 0
    power_loss_of_considered_line = 0
    elictricity_loss = 0
    electricity_loss_cost = 0
    deprecation_value = 0
    capital_investments_switcher = 0
    deprecation_value_for_switcher = 0
    transformer_price = 0
    reduced_losses_in_tp = 0
    loss_electricity_in_tp = 0
    cost_losses_of_electric_energy_in_tp = 0
    deprecation_value_for_tp = 0
    length = 0

    def calc_of_calucated_active_power(self):
        self.calculated_active_power = self.workshop_power * self.demand_coeff

    def calc_of_calucated_reactive_power(self, cosfi):
        self.calculated_reactive_power = self.calculated_active_power * tan(
            acos(cosfi)
        )

    def calc_of_full_power(self):
        self.calculated_full_power = sqrt(
            ((self.calculated_active_power) ** 2)
            + ((self.calculated_reactive_power) ** 2)
        )

    def calc_of_tanfi_n(self):
        self.tan_fi_n = (
            self.calculated_reactive_power / self.calculated_active_power
        )

    def calc_of_compensating_devices(self):
        self.power_of_compensating_devices = self.calculated_active_power * (
            self.tan_fi_n - 0.33
        )

    def compensation_devices_power_choosing(self, transformers_count):
        if (
            self.power_of_compensating_devices <= 25
            and self.transformers_count != 0
        ):
            self.chosen_compens_device = 25
        elif transformers_count != 0:
            self.chosen_compens_device = (
                floor(
                    (
                        (self.power_of_compensating_devices)
                        / (transformers_count * 25)
                    )
                )
                * 25
            )

    def calc_of_complex_reactive_power_after_compensation(
        self, transformers_count
    ):
        self.complex_reactive_power_after_compensation = (
            self.calculated_reactive_power
            - self.transformers_count * self.chosen_compens_device
        )

    def calc_of_full_calculated_transformers_load(self):
        self.full_calculated_transformers_load = sqrt(
            self.calculated_active_power ** 2
            + self.complex_reactive_power_after_compensation ** 2
        )

    def calc_of_calculated_power_of_transformer(self, transformers_count):
        if transformers_count != 0:
            self.calculated_transformers_power = (
                self.full_calculated_transformers_load
                / (transformers_count * 0.7)
            )

    def choosing_transformers_of_workshop(self, transformers_count):
        # if self.full_calculated_transformers_load
        if transformers_count != 0:
            self.preselected_transformer = min(
                n
                for n in nominal_powers_of_transformers
                if n > self.calculated_transformers_power
            )

    def calc_of_transformer_load_coeff(self, transformers_count):
        # for i in range(мощности_трансформаторов):
        if self.transformers_count != 0:
            self.transformer_load_coeff = round(
                (
                    self.calculated_full_power
                    / (
                        self.transformers_count
                        * self.preselected_transformer
                    )
                ),
                2,
            )
            if self.transformer_load_coeff > 0.7:
                self.selected_transformer = nominal_powers_of_transformers[
                    nominal_powers_of_transformers.index(
                        self.preselected_transformer
                    )
                    + 1
                ]
            else:
                self.selected_transformer = (
                    self.preselected_transformer
                )
            # переасчет кз после выбранного трансформатора
            self.transformer_load_coeff = round(
                (
                    self.calculated_full_power
                    / (self.transformers_count * self.selected_transformer)
                ),
                2,
            )
    def calc_of_active_power_loss_ts(self, transformers_types):
        if transformers_types == "масляный" and self.selected_transformer != None:
            self.active_power_loss_ts = delta_p_xx_tm[
                nominal_powers_of_transformers.index(self.selected_transformer)
            ] + delta_p_kz_tm[
                nominal_powers_of_transformers.index(self.selected_transformer)
            ] * pow(
                self.transformer_load_coeff, 2
            )
        elif transformers_types == "сухой" and self.selected_transformer != None:
            self.active_power_loss_ts = delta_p_xx_tsz[
                nominal_powers_of_transformers.index(self.selected_transformer)
            ] + delta_p_kz_tsz[
                nominal_powers_of_transformers.index(self.selected_transformer)
            ] * pow(
                self.transformer_load_coeff, 2
            )
        else:
            pass

    def calc_of_reactive_power_loss_ts(self, transformers_types):
        if transformers_types == "масляный" and self.selected_transformer != None:
            i_xx = i_xx_tm[
                nominal_powers_of_transformers.index(self.selected_transformer)
            ]
            delta_q_xx_tm = (i_xx / 100) * self.selected_transformer
            u_kz = u_kz_tm[
                nominal_powers_of_transformers.index(self.selected_transformer)
            ]
            delta_q_kz_tm = (u_kz / 100) * self.selected_transformer
            self.reactive_power_loss_ts = delta_q_xx_tm + delta_q_kz_tm * pow(
                self.transformer_load_coeff, 2
            )
        elif transformers_types == "сухой" and self.selected_transformer != None:
            i_xx = i_xx_tsz[
                nominal_powers_of_transformers.index(self.selected_transformer)
            ]
            delta_q_xx_tsz = (i_xx / 100) * self.selected_transformer
            u_kz = u_kz_tsz[
                nominal_powers_of_transformers.index(self.selected_transformer)
            ]
            delta_q_kz_tsz = (u_kz / 100) * self.selected_transformer
            self.reactive_power_loss_ts = delta_q_xx_tsz + delta_q_kz_tsz * pow(
                self.transformer_load_coeff, 2
            )

    def расчет_расчетной_нагрузки_а_гпп_тп(self):
        if self.transformers_count != 0:
            self.calculated_load_a_gpp_tp = sqrt(
                (
                    self.calculated_active_power_summ
                    + self.active_power_loss_ts
                )
                ** 2
                + (
                    self.calculated_reactive_power_сумма
                    + self.reactive_power_loss_ts
                )
                ** 2
            )
            print(self.calculated_load_a_gpp_tp, "calculated_load_a_gpp_tp")
            print(
                self.calculated_active_power_summ,
                "calculated_active_power_summ",
            )

    def расчет_расчетной_нагрузки_н_гпп_тп(self):
        if self.transformers_count != 0:
            self.calculated_load_n_gpp_tp = sqrt(
                (
                    (
                        self.calculated_active_power_summ
                        / self.transformers_count
                    )
                    + self.active_power_loss_ts
                )
                ** 2
                + (
                    (
                        self.calculated_reactive_power_сумма
                        / self.transformers_count
                    )
                    + self.reactive_power_loss_ts
                )
                ** 2
            )


    def расчет_расчетной_токовой_нагрузки_а_гпп_тп(self):
        if self.transformers_count != 0:
            self.calculated_current_load_a_gpp_tp = (
                self.calculated_load_a_gpp_tp / (sqrt(3) * 10)
            )
            print(self.calculated_load_a_gpp_tp, "self.calculated_load_a_gpp_tp")
        elif self.transformers_count == 0:
            self.calculated_current_load_a_gpp_tp = (
                self.full_calculated_transformers_load / (sqrt(3) * 10)
            )

    def расчет_расчетной_токовой_нагрузки_н_гпп_тп(self):
        if self.transformers_count != 0:
            self.calculated_current_load_n_gpp_tp = (
                self.calculated_load_n_gpp_tp / (sqrt(3) * 10)
            )
            print(
                self.calculated_current_load_n_gpp_tp,
                "self.calculated_current_load_n_gpp_tp",
            )
        elif self.transformers_count == 0:
            self.calculated_current_load_n_gpp_tp = (
                self.full_calculated_transformers_load / (sqrt(3) * 10)
            )

    def расчет_к1(self, transformers_count):
        if self.transformers_count == 1:
            self.k1 = 1
        elif self.transformers_count == 2:
            self.k1 = 0.9
        elif self.transformers_count == 4:
            self.k1 = 0.8

    def расчет_допустимой_токовой_нагрузка_н_гпп_тп(self):
        if self.transformers_count != 0:
            self.permissible_current_load_n_gpp_tp = (
                self.calculated_current_load_n_gpp_tp / self.k1
            )
        elif self.transformers_count == 0:
            self.permissible_current_load_n_gpp_tp = (
                self.calculated_current_load_n_gpp_tp
            )

    def расчет_допустимой_токовой_нагрузка_а_гпп_тп(self):
        if self.transformers_count != 0:
            self.permissible_current_load_a_gpp_tp = (
                self.calculated_current_load_a_gpp_tp / (self.k1 * self.k2)
            )
        elif self.transformers_count == 0:
            self.permissible_current_load_a_gpp_tp = (
                self.calculated_current_load_a_gpp_tp
            )

    def выбор_кабелей(self, transformers_count):
        if self.transformers_count > 0:
            # cable = None
            try:
                cable = next(
                    num
                    for num in i_dop_kab_aasv
                    if num > self.permissible_current_load_a_gpp_tp
                )
                self.cable_current = cable
                # print(cable, 'cable')

                l_delta_index = i_dop_kab_aasv.index(cable)
                self.permissible_voltage_loss = (
                    l_delta_u1percent[l_delta_index]
                    * 10
                    * (cable / self.permissible_current_load_a_gpp_tp)
                )
                # print(permissible_voltage_loss, 'permissible_voltage_loss')
                while self.permissible_voltage_loss < self.line_lengths:
                    i_dop_kab_aasv_index = i_dop_kab_aasv(cable)
                    cable = i_dop_kab_aasv[i_dop_kab_aasv_index + 1]

                self.cable_current = cable
                # if self.permissible_voltage_loss < self.line_lengths:
                #     i_dop_kab_aasv_index = i_dop_kab_aasv(cable)
                #     cable = i_dop_kab_aasv[i_dop_kab_aasv_index + 1]

                self.permissible_economic_current_density_loss = (
                    self.calculated_load_n_gpp_tp / 1.4
                )
                # print(permissible_economic_current_density_loss, 'permissible_economic_current_density_loss')
                while (
                    self.permissible_economic_current_density_loss
                    > cable_sections[i_dop_kab_aasv.index(cable)]
                ):
                    i_dop_kab_aasv_index = i_dop_kab_aasv.index(cable)
                    cable = i_dop_kab_aasv[i_dop_kab_aasv_index + 1]

                self.cable_current = cable

                self.chosen_cable = cable_sections[i_dop_kab_aasv.index(cable)]
                # print(self.permissible_voltage_loss, 'permissible_voltage_loss')
                # print(self.line_lengths, 'line_lengths')
                # print(self.permissible_current_load_a_gpp_tp, 'permissible_current_load_a_gpp_tp')
                # print(self.permissible_economic_current_density_loss, 'permissible_economic_current_density_loss')
            except StopIteration:
                print("done")
        elif self.transformers_count == 0:
            try:
                cable = next(
                    num
                    for num in i_dop_kab_asb
                    if num > self.permissible_current_load_a_gpp_tp
                )
                self.cable_current = cable
                i_dop_kab_asb_index = i_dop_kab_asb.index(cable)
                if (
                    self.full_calculated_transformers_load
                    < s_asb[i_dop_kab_asb_index]
                ):
                    cable = i_dop_kab_asb[i_dop_kab_asb_index]
                    self.cable_current = cable
                else:
                    i_dop_kab_asb_index = i_dop_kab_asb.index(cable)
                    while (
                        self.full_calculated_transformers_load
                        < s_asb[i_dop_kab_asb_index]
                    ):
                        i_dop_kab_asb_index = i_dop_kab_asb.index(cable)
                        cable = i_dop_kab_asb[i_dop_kab_asb_index + 1]
                    self.cable_current = cable

                self.cable_current = cable

                self.chosen_cable = cable_sections[i_dop_kab_asb.index(cable)]

                print(self.chosen_cable, "cable")

            except StopIteration:
                print("done")

    def расчет_капитальных_вложений_на_линию(
        self, transformers_count, line_lengths
    ):
        if transformers_count != 0:
            self.line_investments = (
                transformers_count
                * line_lengths
                * pow(10, -3)
                * one_km_cost_asb_cost_one_km_line_aasv[
                    cable_sections.index(self.chosen_cable)
                ]
            )
        elif transformers_count == 0:
            self.line_investments = (
                line_lengths
                * pow(10, -3)
                * one_km_cost_asb_cost_one_km_line_asb[
                    cable_sections.index(self.chosen_cable)
                ]
            )

    def расчет_расхода_цветного_материала(
        self, transformers_count, line_lengths
    ):
        if transformers_count > 0:
            self.color_material_consumption = (
                self.transformers_count
                * line_lengths
                * pow(10, -3)
                * specific_consumption_of_color_material_aashv[
                    cable_sections.index(self.chosen_cable)
                ]
            )
        elif transformers_count == 0:
            self.color_material_consumption = (
                line_lengths
                * pow(10, -3)
                * specific_color_material_consumption_asb[
                    cable_sections.index(self.chosen_cable)
                ]
            )

    def расчет_коэффициента_загрузки(self, transformers_count):
        # if transformers_count > 0:
        self.load_coeff = (
            self.calculated_current_load_n_gpp_tp / self.cable_current
        )

    def расчет_потерь_мощности_рассматриваемой_линии(
        self, transformers_count, line_lengths
    ):
        if transformers_count > 0:
            self.power_loss_of_considered_line = (
                transformers_count
                * self.active_power_loss_ts
                * line_lengths
                * pow(10, -3)
                * self.load_coeff ** 2
            )
        # elif transformers_count == 0:
        #     self.power_loss_of_considered_line = self.active_power_loss_ts * line_lengths * pow(10, -3) * self.load_coeff**2

    def расчет_потерь_электроэнергии(self):
        self.elictricity_loss = (
            self.power_loss_of_considered_line * count_of_working_hours
        )

    def расчет_стоимости_потерь_электроэнергии(self):
        self.electricity_loss_cost = (
            self.elictricity_loss * cost__of_kw_electricity * pow(10, -3)
        )

    def расчет_стоимости_амортизационных_отчислений(self):
        self.deprecation_value = (
            coefficient_depreciation_deductions * self.line_investments
        )

    def расчет_капитальных_вложений_выключатель(
        self, switchers_count, switcher_prices
    ):
        self.capital_investments_switcher = (
            switchers_count * switcher_prices
        )

    def расчет_стоимости_амортизационных_отчислений_выключатель(self):
        self.deprecation_value_for  = (
            coefficient_depreciation_deductions_switcher
            * self.capital_investments_switcher
        )

    def расчет_стоимости_трансформаторов(
        self, transformers_count
    ):
        if transformers_count != 0:
            self.transformer_price = (
                transformers_count * transformers_prices[nominal_powers_of_transformers.index(self.selected_transformer)]
            )

    def расчет_привиденных_потерь_в_тп(
        self, transformers_count, transformers_types
    ):
        if transformers_types == "масляный" and transformers_count != 0:
            self.reduced_losses_in_tp = transformers_count * (
                delta_p_xx_tm[
                    nominal_powers_of_transformers.index(
                        self.selected_transformer
                    )
                ]
                + delta_p_kz_tm[
                    nominal_powers_of_transformers.index(
                        self.selected_transformer
                    )
                ]
                * self.transformer_load_coeff ** 2
            )
        elif transformers_types == "сухой" and transformers_count != 0:
            self.reduced_losses_in_tp = transformers_count * (
                delta_p_xx_tsz[
                    nominal_powers_of_transformers.index(
                        self.selected_transformer
                    )
                ]
                + delta_p_kz_tsz[
                    nominal_powers_of_transformers.index(
                        self.selected_transformer
                    )
                ]
                * self.transformer_load_coeff ** 2
            )

    def расчет_потерь_электроэнергии_в_тп(
        self, transformers_count, transformers_types
    ):
        if transformers_types == "масляный" and transformers_count != 0:
            self.loss_electricity_in_tp = transformers_count * (
                delta_p_xx_tm[
                    nominal_powers_of_transformers.index(
                        self.selected_transformer
                    )
                ]
                * т_в
                + delta_p_kz_tm[
                    nominal_powers_of_transformers.index(
                        self.selected_transformer
                    )
                ]
                * (self.transformer_load_coeff ** 2)
                * count_of_working_hours
            )
        elif transformers_types == "сухой" and transformers_count != 0:
            self.loss_electricity_in_tp = transformers_count * (
                delta_p_xx_tsz[
                    nominal_powers_of_transformers.index(
                        self.selected_transformer
                    )
                ]
                * т_в
                + delta_p_kz_tsz[
                    nominal_powers_of_transformers.index(
                        self.selected_transformer
                    )
                ]
                * (self.transformer_load_coeff ** 2)
                * count_of_working_hours
            )

    def расчет_стоимости_потерь_электроэнергии_в_тп(self, transformers_count):
        if transformers_count != 0:
            self.cost_losses_of_electric_energy_in_tp = (
                self.loss_electricity_in_tp
                * cost__of_kw_electricity
                * pow(10, -3)
            )

    def расчет_стоимости_амортизационных_отчислений_на_тп(
        self, transformers_count
    ):
        if transformers_count != 0:
            self.deprecation_value_for_tp = (
                transformers_count
                * transformers_prices[nominal_powers_of_transformers.index(self.selected_transformer)]
                * coefficient_depreciation_deductions_tp
            )

    def расчет_длины(self, workshops_coords):
        self.length = abs(gpp_coords[0] - workshops_coords[0]) + abs(
            gpp_coords[1] - abs(workshops_coords[1])
        )

    def __init__(
        self,
        workshop_power,
        demand_coeff,
        cos_fi,
        transformers_count,
        transformers_types,
        line_lengths,
        switchers_count,
        switcher_types,
        switcher_prices,
        workshops_coords,
    ):
        self.workshop_power = workshop_power
        self.demand_coeff = demand_coeff
        self.cos_fi = cos_fi
        self.transformers_count = transformers_count
        self.transformers_types = transformers_types
        self.line_lengths = line_lengths
        self.switchers_count = switchers_count
        self.switcher_prices = switcher_prices
        self.switcher_types = switcher_types
        self.workshops_coords = workshops_coords

        self.calc_of_calucated_active_power()
        self.calc_of_calucated_reactive_power(cos_fi)
        self.calc_of_full_power()
        self.calc_of_tanfi_n()
        self.calc_of_compensating_devices()
        self.compensation_devices_power_choosing(transformers_count)
        self.calc_of_complex_reactive_power_after_compensation(
            transformers_count
        )
        self.calc_of_full_calculated_transformers_load()
        self.calc_of_calculated_power_of_transformer(transformers_count)
        self.choosing_transformers_of_workshop(transformers_count)
        self.calc_of_transformer_load_coeff(transformers_count)
        self.calc_of_active_power_loss_ts(transformers_types)
        self.calc_of_reactive_power_loss_ts(transformers_types)
        self.расчет_расчетной_нагрузки_а_гпп_тп()
        self.расчет_расчетной_нагрузки_н_гпп_тп()
        self.расчет_к1(transformers_count)
        self.расчет_расчетной_токовой_нагрузки_а_гпп_тп()
        self.расчет_расчетной_токовой_нагрузки_н_гпп_тп()
        self.расчет_допустимой_токовой_нагрузка_н_гпп_тп()
        self.расчет_допустимой_токовой_нагрузка_а_гпп_тп()
        self.выбор_кабелей(transformers_count)
        self.расчет_капитальных_вложений_на_линию(
            transformers_count, line_lengths
        )
        self.расчет_расхода_цветного_материала(transformers_count, line_lengths)
        self.расчет_коэффициента_загрузки(transformers_count)
        self.расчет_потерь_мощности_рассматриваемой_линии(
            transformers_count, line_lengths
        )
        self.расчет_потерь_электроэнергии()
        self.расчет_стоимости_потерь_электроэнергии()
        self.расчет_стоимости_амортизационных_отчислений()
        self.расчет_капитальных_вложений_выключатель(
            switchers_count, switcher_prices
        )
        self.расчет_стоимости_трансформаторов(
            transformers_count
        )
        self.расчет_привиденных_потерь_в_тп(
            transformers_count, transformers_types
        )
        self.расчет_потерь_электроэнергии_в_тп(
            transformers_count, transformers_types
        )
        self.расчет_стоимости_потерь_электроэнергии_в_тп(transformers_count)
        self.расчет_стоимости_амортизационных_отчислений_на_тп(
            transformers_count
        )
        self.расчет_длины(workshops_coords)

    # def __repr__(self):
    # return "<Цех input_param:%d some_param:%d other_param:%d>" %
    # (self.input_param, self.some_param, self.other_param)

    # def __str__(self):
    # return "<Цех input_param:%d some_param:%d other_param:%d>" %
    # (self.input_param, self.some_param, self.other_param)


нет_доп_цехов = -1
карта_дополнительных_цехов = {
    0: [нет_доп_цехов],
    1: [нет_доп_цехов],
    2: [нет_доп_цехов],
    3: [нет_доп_цехов],
    4: [нет_доп_цехов],
    5: [нет_доп_цехов],
    6: [нет_доп_цехов],
    7: [нет_доп_цехов],
    8: [6],
    9: [нет_доп_цехов],
    10: [нет_доп_цехов],
    11: [нет_доп_цехов],
    12: [11],
    13: [10],
    14: [нет_доп_цехов],
}  # НУЖНО НЕ ЗАБЫТЬ УЧЕСТЬ, что это индексы цехов, а не сами цеха


# noinspection NonAsciiCharacters
class КонтейнерЦехов:
    line_length = 0.6
    max_time_difference_coeff = 0.9
    активная_мощность_с_учетом_потерь = 0
    реактивная_мощность_с_учетом_потерь = 0
    power_of_compensating_devices = 0
    потери_мощности_в_ку = 0
    calculated_active_power_on_gpp = 0
    calculated_reactive_power_on_gpp = 0
    полная_расчетная_мощность_на_шинах_гпп = 0
    потери_активной_мощности_в_трансф_гпп = 0
    потери_реактивной_мощности_в_трансф_гпп = 0
    полная_расчетная_нагрузка_с_учетом_потерь_мощности = 0
    напряжение_питающей_линии = 0

    # к_л = 0
    # к_эа = 0
    # к_тп = 0
    # с_ал = 0
    # с_пт = 0
    # с_пл = 0
    # с_аэа = 0
    # с_атп = 0
    # затраты = 0

    def расчет_активной_мощности_с_учетом_потерь(self):
        сумма_расчётных_активных_мощностей = 0

        for цех in self.цеха:
            сумма_расчётных_активных_мощностей += цех.calculated_active_power
        self.активная_мощность_с_учетом_потерь = сумма_расчётных_активных_мощностей

    def расчет_реактивной_мощности_с_учетом_потерь(self):
        сумма_расчётных_реактивных_мощностей = 0

        for цех in self.цеха:
            сумма_расчётных_реактивных_мощностей += цех.calculated_reactive_power
        self.реактивная_мощность_с_учетом_потерь = сумма_расчётных_реактивных_мощностей

    def calc_of_compensating_devices(self):
        tanfi_n = (
            self.реактивная_мощность_с_учетом_потерь
            / self.активная_мощность_с_учетом_потерь
        )
        # print(tanfi_n, 'tanfi_n')
        self.power_of_compensating_devices = self.активная_мощность_с_учетом_потерь * (
            tanfi_n - 0.33
        )

    def расчет_потери_мощности_в_ку(self):
        self.потери_мощности_в_ку = 0.002 * self.power_of_compensating_devices

    def расчет_расчетной_активной_мощности_на_шинах_гпп(self):
        self.calculated_active_power_on_gpp = (
            self.активная_мощность_с_учетом_потерь + self.потери_мощности_в_ку
        ) * self.max_time_difference_coeff

    def расчет_расчетной_реактивной_мощности_на_шинах_гпп(self):
        self.calculated_reactive_power_on_gpp = (
            self.реактивная_мощность_с_учетом_потерь - self.power_of_compensating_devices
        ) * self.max_time_difference_coeff

    def расчет_полной_расчетной_мощности_на_шинах_гпп(self):
        self.полная_расчетная_мощность_на_шинах_гпп = sqrt(
            self.calculated_active_power_on_gpp ** 2
            + self.calculated_reactive_power_on_gpp ** 2
        )

    def расчет_потерь_активной_мощности_в_трансф_гпп(self):
        self.потери_активной_мощности_в_трансф_гпп = (
            0.02 * self.полная_расчетная_мощность_на_шинах_гпп
        )

    def расчет_потерь_реактивной_мощности_в_трансф_гпп(self):
        self.потери_реактивной_мощности_в_трансф_гпп = (
            0.1 * self.полная_расчетная_мощность_на_шинах_гпп
        )

    def расчет_полной_расчетной_нагрузки_с_учетом_потерь_мощности(self):
        self.полная_расчетная_нагрузка_с_учетом_потерь_мощности = sqrt(
            (
                self.calculated_active_power_on_gpp
                + self.потери_активной_мощности_в_трансф_гпп
            )
            ** 2
            + (
                self.calculated_reactive_power_on_gpp
                + self.потери_реактивной_мощности_в_трансф_гпп
            )
            ** 2
        )

    def расчет_напряжения_питающей_линии(self):
        u1 = (
            3 * sqrt(self.полная_расчетная_нагрузка_с_учетом_потерь_мощности * 10 ** -3)
            + 0.5 * self.line_length
        )

        u2 = 4.34 * sqrt(
            self.line_length
            + 16
            * (
                self.calculated_active_power_on_gpp
                + self.потери_активной_мощности_в_трансф_гпп
            )
            * 10 ** -3
        )
        u3 = 16 * sqrt(
            sqrt(
                (
                    self.calculated_active_power_on_gpp
                    + self.потери_активной_мощности_в_трансф_гпп
                )
                * 10 ** -3
                * self.line_length
            )
        )
        u4 = 17 * sqrt(
            (self.line_length / 16)
            + (
                self.calculated_active_power_on_gpp
                + self.потери_активной_мощности_в_трансф_гпп
            )
            * 10 ** -3
        )

        sum_u = u1 + u2 + u3 + u4
        self.напряжение_питающей_линии = sum_u / 4

    def расчёт_дополнительной_мощности(
        self, основной_цех: Цех, дополнительные_цеха: [Цех]
    ):
        основной_цех.calculated_reactive_power_сумма = (
            основной_цех.calculated_reactive_power
        )
        основной_цех.calculated_active_power_summ = (
            основной_цех.calculated_active_power
        )

        for дополнительный_цех in дополнительные_цеха:
            основной_цех.calculated_reactive_power_сумма += (
                дополнительный_цех.calculated_reactive_power
            )
            основной_цех.calculated_active_power_summ += (
                дополнительный_цех.calculated_active_power
            )

    к_л = 0
    к_эа = 0
    к_тп = 0
    с_ал = 0
    с_пт = 0
    с_пл = 0
    с_аэа = 0
    с_атп = 0
    затраты = 0

    def расчет__затрат(self):
        for i in range(workshops_counts):
            self.к_л += self.цеха[i].line_investments
            self.к_эа += self.цеха[i].capital_investments_switcher
            self.к_тп += self.цеха[i].transformer_price
            self.с_ал += self.цеха[i].deprecation_value
            self.с_пт += self.цеха[i].cost_losses_of_electric_energy_in_tp
            self.с_пл += self.цеха[i].electricity_loss_cost
            self.с_аэа += self.цеха[i].deprecation_value_for_switcher
            self.с_атп += self.цеха[i].deprecation_value_for_tp

    def расчет_затрат(self):
        self.затраты = (
            0.15 * (self.к_л + self.к_эа + self.к_тп)
            + self.с_ал
            + self.с_пт
            + self.с_пл
            + self.с_аэа
            + self.с_атп
        )

    def __init__(self):
        self.цеха = []
        for i in range(len(enterprise_powers)):
            цех = Цех(
                enterprise_powers[i],
                demand_coefficients[i],
                cos_fi[i],
                transformers_count[i],
                # мощности_трансформаторов[i],
                transformers_types[i],
                line_lengths[i],
                switcher_types[i],
                switcher_prices[i],
                switchers_count[i],
                workshops_coords[i],
            )
            self.цеха.append(цех)

        # print(self.цеха[0].calculated_active_power  , 'calculated_active_power')
        self.расчет_активной_мощности_с_учетом_потерь()
        self.расчет_реактивной_мощности_с_учетом_потерь()
        # print(self.активная_мощность_с_учетом_потерь, 'активная_мощность_с_учетом_потерь')
        # print(self.реактивная_мощность_с_учетом_потерь, 'реактивная_мощность_с_учетом_потерь')
        self.calc_of_compensating_devices()
        self.расчет_потери_мощности_в_ку()
        self.расчет_расчетной_активной_мощности_на_шинах_гпп()
        self.расчет_расчетной_реактивной_мощности_на_шинах_гпп()
        self.расчет_полной_расчетной_мощности_на_шинах_гпп()
        self.расчет_потерь_активной_мощности_в_трансф_гпп()
        self.расчет_потерь_реактивной_мощности_в_трансф_гпп()
        self.расчет_полной_расчетной_нагрузки_с_учетом_потерь_мощности()
        self.расчет_напряжения_питающей_линии()
        self.расчет__затрат()
        self.расчет_затрат()

        # учет РП(мощности тп + мощности тп)

        for _, (индекс_основного_цех, индексы_дополнительных_цехов) in enumerate(
            карта_дополнительных_цехов.items()
        ):
            основной_цех = self.цеха[индекс_основного_цех]
            дополнительные_цеха = []

            if индексы_дополнительных_цехов[0] == нет_доп_цехов:
                основной_цех.calculated_reactive_power_сумма = (
                    основной_цех.calculated_reactive_power
                )
                основной_цех.calculated_active_power_summ = (
                    основной_цех.calculated_active_power
                )
                continue

            for индекс_дополнительного_цеха in индексы_дополнительных_цехов:
                дополнительные_цеха.append(self.цеха[индекс_дополнительного_цеха])

            self.расчёт_дополнительной_мощности(
                основной_цех=основной_цех, дополнительные_цеха=дополнительные_цеха
            )


# noinspection NonAsciiCharacters
def main():
    workshops_conatiner = КонтейнерЦехов()

    рп1 = "1 рп"
    рп2 = "2 рп"
    тр1 = "1 тр"
    тр2 = "2 тр"

    категории = [2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 3, 3, 2, 2, 2]

    варианты_схем = [
        [рп1, тр1, тр2, рп2] if категория_цехов == 3 else [тр2]
        for категория_цехов in категории
    ]

    for index, категория_цехов in enumerate(workshops_categories):
        if категория_цехов == 3:
            if (
                workshops_conatiner.цеха[index].full_calculated_transformers_load
                <= 330
            ):
                pass

            elif (
                workshops_conatiner.цеха[index].full_calculated_transformers_load
                > 330
                and workshops_conatiner.цеха[
                    index
                ].full_calculated_transformers_load
                <= 660
            ):

                варианты_схем[index].remove(рп1)

            elif (
                workshops_conatiner.цеха[index].full_calculated_transformers_load
                > 660
            ):
                варианты_схем[index].remove(рп1)

                варианты_схем[index].remove(рп2)

    # print(варианты_схем)
    # ------------------------s
    список_всех_возможных_варианций_схем_с_учетом = list(
        itertools.product(*варианты_схем)
    )
    workshops_counts = 15
    # ------------------------e

    # ------------------------s
    # arr = []

    # for вариант_схемы in список_всех_возможных_варианций_схем_с_учетом:
    #     # print(вариант_схемы, 'вариант_схемы')
    #     нули_и_единицы = [0] * workshops_counts
    #     for цех in range(len(вариант_схемы)):
    #         if вариант_схемы[цех] == тр2 or вариант_схемы[цех] == тр1:
    #             нули_и_единицы[цех] = 1
    #         else:
    #             нули_и_единицы[цех] = 0
    #     arr.append(нули_и_единицы)
    # ------------------------e

    # print(arr[2], 'arr[2]')
    # array.append(zeros_and_ones)
    # print(arr[2])

    # ------------------------s
    # for вариант_схемы in arr:
    #     for цех in range(len(вариант_схемы)):
    #         if вариант_схемы[цех] == 1:
    #             вариант_схемы[цех] = [-1]
    #         elif вариант_схемы[цех] == 0:
    #             вариант_схемы[цех] = []

    #     for цех in range(len(вариант_схемы)):
    #         if вариант_схемы[цех] == []:
    #             for цех_изб in range(len(вариант_схемы)):
    #                 if вариант_схемы[цех_изб] == [-1]:
    #                     вариант_схемы[цех].append(цех_изб)
    # ------------------------e

    # print(arr[200], 'arr 200')
    # print(список_всех_возможных_варианций_схем_с_учетом[500], 'индекс 500')
    # print(arr[500], 'индекс 500')
    # print(arr[20], 'индекс 20')
    # print(список_всех_возможных_варианций_схем_с_учетом[20], 'индекс 20')
    # print(arr[2], 'index 1')
    # print(arr[5000], 'index 2')
    # print(zeros_and_ones)
    # ar = [[[-1], [-1], [-1], [-1], [-1], [-1], [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 14], [-1], [-1], [-1], [-1], [-1], [-1], [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 14], [-1], [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 14]],
    # [[-1], [-1], [-1], [-1], [-1], [-1], [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11], [-1], [-1], [-1], [-1], [-1], [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11], [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11], [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11], [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11]]]

    # ------------------------s
    # список_всех_возможных_подключений = []
    # for список in arr:
    #     список_всех_возможных_подключений.append(
    #         list(itertools.product(*список)))
    # print(список_всех_возможных_подключений)
    # ------------------------e

    # список_всех_возможных_подключений = list(itertools.product(*ar[0]))
    # print(список_всех_возможных_подключений[1000])
    # print(список_всех_возможных_подключений)
    # result = []
    # for p in список_всех_возможных_подключений:
    #     result.append(p)
    # print(result)
    # print(len(result))
    # # for key, value in карта_доп_цехов:

    # варианты_подключений = []

    # for index, категория_цехов in enumerate(workshops_categories):
    #     if категория_цехов == 3:
    #         for i in range(len(workshops_categories)):
    #             вариант_подключения = {}
    #             for key, val in вариант_подключения:
    #                 вариант_подключения[key] = [index]

    # варианты_подключений[i] =

    # print(workshops_conatiner.потери_мощности_в_ку)
    # print(workshops_conatiner.полная_расчетная_мощность_на_шинах_гпп)
    # print(workshops_conatiner.полная_расчетная_нагрузка_с_учетом_потерь_мощности)
    # print(workshops_conatiner.напряжение_питающей_линии, 'напряжение_питающей_линии')
    # print(workshops_conatiner.цеха[0].calculated_active_power, 'calculated_active_power цех 1')
    # print(workshops_conatiner.цеха[0].calculated_reactive_power, 'calculated_reactive_power цех 1')
    # print(workshops_conatiner.цеха[0].calculated_active_power_summ, 'calculated_active_power_summ цех 1')
    # print(workshops_conatiner.цеха[0].calculated_reactive_power_сумма, 'calculated_reactive_power_сумма цех 1')
    # print(workshops_conatiner.цеха[3].calculated_active_power, 'calculated_active_power цех 4')
    # print(workshops_conatiner.цеха[6].calculated_active_power, 'calculated_active_power цех 7 ')
    # print(workshops_conatiner.цеха[9].calculated_active_power, 'calculated_active_power цех 10 ')
    # print(workshops_conatiner.цеха[3].calculated_active_power_summ, 'calculated_active_power_summ')
    # print(workshops_conatiner.цеха[0].transformer_load_coeff, 'transformer_load_coeff')
    # print(workshops_conatiner.цеха[0].reactive_power_loss_ts,
    #       'reactive_power_loss_ts')
    # print(workshops_conatiner.цеха[0].calculated_load_a_gpp_tp,
    #       'calculated_load_a_gpp_tp')
    # print(workshops_conatiner.цеха[0].calculated_load_n_gpp_tp,
    #       'calculated_load_n_gpp_tp')
    # print(workshops_conatiner.цеха[1].calculated_load_a_gpp_tp, 'calculated_load_a_gpp_tp')
    # print(workshops_conatiner.цеха[1].calculated_load_n_gpp_tp, 'calculated_load_n_gpp_tp')
    # print(workshops_conatiner.цеха[0].calculated_current_load_n_gpp_tp, 'calculated_current_load_n_gpp_tp')
    # print(workshops_conatiner.цеха[1].permissible_voltage_loss, 'permissible_voltage_loss')
    # print(workshops_conatiner.цеха[1].permissible_economic_current_density_loss, 'permissible_economic_current_density_loss')
    # print(workshops_conatiner.цеха[0].permissible_current_load_a_gpp_tp, 'permissible_current_load_a_gpp_tp')
    # print(workshops_conatiner.цеха[6].full_calculated_transformers_load, 'full_calculated_transformers_load')
    # print(workshops_conatiner.цеха[3].calculated_active_power_summ, 'calculated_active_power_summ')
    # print(workshops_conatiner.цеха[3].calculated_active_power, 'calculated_active_power')
    # print(workshops_conatiner.цеха[6].calculated_active_power, 'calculated_active_power')
    # print(workshops_conatiner.цеха[0].transformer_load_coeff, 'transformer_load_coeff')
    print(
        workshops_conatiner.цеха[0].electricity_loss_cost,
        "electricity_loss_cost",
    )
    # print(workshops_conatiner.цеха[0].cost_losses_of_electric_energy_in_tp, 'cost_losses_of_electric_energy_in_tp')
    print(
        workshops_conatiner.цеха[11].transformer_load_coeff,
        "transformer_load_coeff",
    )
    print(
        workshops_conatiner.цеха[11].full_calculated_transformers_load,
        "full_calculated_transformers_load",
    )
    print(
        workshops_conatiner.цеха[11].calculated_transformers_power,
        "calculated_transformers_power",
    )
    print(
        workshops_conatiner.цеха[11].preselected_transformer,
        "предварительно_selected_transformer",
    )
    print(workshops_conatiner.цеха[11].selected_transformer, "selected_transformer")
    print(workshops_conatiner.затраты, "затраты")
    # print(workshops_conatiner.цеха[9].length, 'length')
    # print(workshops_conatiner.цеха[1].length, 'length')
    # print(workshops_conatiner.цеха[0].complex_reactive_power_after_compensation, 'complex_reactive_power_after_compensation')

    # print(workshops_conatiner.цеха[1].chosen_cable, 'chosen_cable ЦЕХ 2')
    # print(workshops_conatiner.цеха[0].color_material_consumption, 'color_material_consumption')
    # print(workshops_conatiner.цеха[0].cable_current, 'cable_current')

    # print(workshops_conatiner.цеха[0].load_coeff, 'коэффициент_загрузк��')
    # for цех in workshops_conatiner.цеха:
    #     print(цех.full_calculated_transformers_load, 'full_calculated_transformers_load ')


if __name__ == "__main__":
    main()
