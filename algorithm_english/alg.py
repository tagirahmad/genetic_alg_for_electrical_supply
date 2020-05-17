# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-

from math import sqrt, ceil, acos, tan, floor
import itertools
from constants.constants import *
from constants.aasv_cable import *
from constants.asb_cable import *
from constants.circuit_breaker import *
from constants.transformer import *


class Workshop:
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

    def calc_of_calculated_load_a_gpp_tp(self):
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

    def calc_of_calculated_load_n_gpp_tp(self):
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

    def calc_of_k1(self, transformers_count):
        if self.transformers_count == 1:
            self.k1 = 1
        elif self.transformers_count == 2:
            self.k1 = 0.9
        elif self.transformers_count == 4:
            self.k1 = 0.8

    def calc_of_permissible_current_load_n_gpp_tp(self):
        if self.transformers_count != 0:
            self.permissible_current_load_n_gpp_tp = (
                self.calculated_current_load_n_gpp_tp / self.k1
            )
        elif self.transformers_count == 0:
            self.permissible_current_load_n_gpp_tp = (
                self.calculated_current_load_n_gpp_tp
            )

    def calc_of_permissible_current_load_a_gpp_tp(self):
        if self.transformers_count != 0:
            self.permissible_current_load_a_gpp_tp = (
                self.calculated_current_load_a_gpp_tp / (self.k1 * self.k2)
            )
        elif self.transformers_count == 0:
            self.permissible_current_load_a_gpp_tp = (
                self.calculated_current_load_a_gpp_tp
            )

    def cable_choosing(self, transformers_count):
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

    def calc_of_line_investments(
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

    def calc_of_color_material_consumption(
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

    def calc_of_load_coeff(self, transformers_count):
        if transformers_count > 0:
            self.load_coeff = (
                self.calculated_current_load_n_gpp_tp / self.cable_current
            )

    def calc_of_power_loss_of_considered_line(
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

    def calc_of_elictricity_loss(self):
        self.elictricity_loss = (
            self.power_loss_of_considered_line * count_of_working_hours
        )

    def calc_of_electricity_loss_cost(self):
        self.electricity_loss_cost = (
            self.elictricity_loss * cost__of_kw_electricity * pow(10, -3)
        )

    def calc_of_deprecation_value(self):
        self.deprecation_value = (
            coefficient_depreciation_deductions * self.line_investments
        )

    def calc_of_capital_investments_switcher(
        self, switchers_count, switcher_prices
    ):
        self.capital_investments_switcher = (
            switchers_count * switcher_prices
        )

    def calc_of_deprecation_value_for_switcher(self):
        self.deprecation_value_for_switcher  = (
            coefficient_depreciation_deductions_switcher
            * self.capital_investments_switcher
        )

    def calc_of_transformer_price(
        self, transformers_count
    ):
        if transformers_count != 0:
            self.transformer_price = (
                transformers_count * transformers_prices[nominal_powers_of_transformers.index(self.selected_transformer)]
            )

    def calc_of_reduced_losses_in_tp(
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

    def calc_of_loss_electricity_in_tp(
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

    def calc_of_cost_losses_of_electric_energy_in_tp(self, transformers_count):
        if transformers_count != 0:
            self.cost_losses_of_electric_energy_in_tp = (
                self.loss_electricity_in_tp
                * cost__of_kw_electricity
                * pow(10, -3)
            )

    def calc_of_deprecation_value_for_tp(
        self, transformers_count
    ):
        if transformers_count != 0:
            self.deprecation_value_for_tp = (
                transformers_count
                * transformers_prices[nominal_powers_of_transformers.index(self.selected_transformer)]
                * coefficient_depreciation_deductions_tp
            )

    def calc_of_length(self, workshops_coords):
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
        self.calc_of_calculated_load_a_gpp_tp()
        self.calc_of_calculated_load_n_gpp_tp()
        self.calc_of_k1(transformers_count)
        self.расчет_расчетной_токовой_нагрузки_а_гпп_тп()
        self.расчет_расчетной_токовой_нагрузки_н_гпп_тп()
        self.calc_of_permissible_current_load_n_gpp_tp()
        self.calc_of_permissible_current_load_a_gpp_tp()
        self.cable_choosing(transformers_count)
        self.calc_of_line_investments(
            transformers_count, line_lengths
        )
        self.calc_of_color_material_consumption(transformers_count, line_lengths)
        self.calc_of_load_coeff(transformers_count)
        self.calc_of_power_loss_of_considered_line(
            transformers_count, line_lengths
        )
        self.calc_of_elictricity_loss()
        self.calc_of_electricity_loss_cost()
        self.calc_of_deprecation_value()
        self.calc_of_capital_investments_switcher(
            switchers_count, switcher_prices
        )
        self.calc_of_deprecation_value_for_switcher()
        self.calc_of_transformer_price(
            transformers_count
        )
        self.calc_of_reduced_losses_in_tp(
            transformers_count, transformers_types
        )
        self.calc_of_loss_electricity_in_tp(
            transformers_count, transformers_types
        )
        self.calc_of_cost_losses_of_electric_energy_in_tp(transformers_count)
        self.calc_of_deprecation_value_for_tp(
            transformers_count
        )
        self.calc_of_length(workshops_coords)

    # def __repr__(self):
    # return "<workshop input_param:%d some_param:%d other_param:%d>" %
    # (self.input_param, self.some_param, self.other_param)

    # def __str__(self):
    # return "<workshop input_param:%d some_param:%d other_param:%d>" %
    # (self.input_param, self.some_param, self.other_param)


no_additional_workshops = -1
additional_workshops_map = {
    0: [no_additional_workshops],
    1: [no_additional_workshops],
    2: [no_additional_workshops],
    3: [no_additional_workshops],
    4: [no_additional_workshops],
    5: [no_additional_workshops],
    6: [no_additional_workshops],
    7: [no_additional_workshops],
    8: [6],
    9: [no_additional_workshops],
    10: [no_additional_workshops],
    11: [no_additional_workshops],
    12: [11],
    13: [10],
    14: [no_additional_workshops],
}  # НУЖНО НЕ ЗАБЫТЬ УЧЕСТЬ, что это индексы цехов, а не сами цеха


# noinspection NonAsciiCharacters
class WorkshopsContainer:
    line_length = 0.6
    max_time_difference_coeff = 0.9
    active_power_with_loss_account = 0
    reactive_power_with_loss_account = 0
    power_of_compensating_devices = 0
    power_loss_in_compensational_devices = 0
    calculated_active_power_on_gpp = 0
    calculated_reactive_power_on_gpp = 0
    full_calculated_power_on_tires_gpp = 0
    active_power_loss_in_transf_gpp = 0
    reactive_power_loss_in_transf_gpp = 0
    full_calculated_load_with_account_of_power_losses = 0
    supply_line_voltage = 0

    # к_л = 0
    # к_эа = 0
    # к_тп = 0
    # с_ал = 0
    # с_пт = 0
    # с_пл = 0
    # с_аэа = 0
    # с_атп = 0
    # expenses = 0

    def calc_of_active_power_with_loss_account(self):
        calculated_active_power_sum = 0

        for workshop in self.workshops:
            calculated_active_power_sum += workshop.calculated_active_power
        self.active_power_with_loss_account = calculated_active_power_sum

    def calc_of_reactive_power_with_loss_account(self):
        calculated_reactive_power_sum = 0

        for workshop in self.workshops:
            calculated_reactive_power_sum += workshop.calculated_reactive_power
        self.reactive_power_with_loss_account = calculated_reactive_power_sum

    def calc_of_compensating_devices(self):
        tanfi_n = (
            self.reactive_power_with_loss_account
            / self.active_power_with_loss_account
        )
        # print(tanfi_n, 'tanfi_n')
        self.power_of_compensating_devices = self.active_power_with_loss_account * (
            tanfi_n - 0.33
        )

    def calc_of_power_loss_in_compensational_devices(self):
        self.power_loss_in_compensational_devices = 0.002 * self.power_of_compensating_devices

    def calc_of_calculated_active_power_on_gpp(self):
        self.calculated_active_power_on_gpp = (
            self.active_power_with_loss_account + self.power_loss_in_compensational_devices
        ) * self.max_time_difference_coeff

    def calc_of_calculated_reactive_power_on_tires_gpp(self):
        self.calculated_reactive_power_on_gpp = (
            self.reactive_power_with_loss_account - self.power_of_compensating_devices
        ) * self.max_time_difference_coeff

    def calc_of_calculated_full_power_on_tires_gpp(self):
        self.full_calculated_power_on_tires_gpp = sqrt(
            self.calculated_active_power_on_gpp ** 2
            + self.calculated_reactive_power_on_gpp ** 2
        )

    def calc_of_loss_active_power_in_transf_gpp(self):
        self.active_power_loss_in_transf_gpp = (
            0.02 * self.full_calculated_power_on_tires_gpp
        )

    def calc_of_loss_reactive_power_in_transf_gpp(self):
        self.reactive_power_loss_in_transf_gpp = (
            0.1 * self.full_calculated_power_on_tires_gpp
        )

    def calc_of_full_calculated_load_with_account_of_power_losses(self):
        self.full_calculated_load_with_account_of_power_losses = sqrt(
            (
                self.calculated_active_power_on_gpp
                + self.active_power_loss_in_transf_gpp
            )
            ** 2
            + (
                self.calculated_reactive_power_on_gpp
                + self.reactive_power_loss_in_transf_gpp
            )
            ** 2
        )

    def calc_of_supply_line_voltage(self):
        u1 = (
            3 * sqrt(self.full_calculated_load_with_account_of_power_losses * 10 ** -3)
            + 0.5 * self.line_length
        )

        u2 = 4.34 * sqrt(
            self.line_length
            + 16
            * (
                self.calculated_active_power_on_gpp
                + self.active_power_loss_in_transf_gpp
            )
            * 10 ** -3
        )
        u3 = 16 * sqrt(
            sqrt(
                (
                    self.calculated_active_power_on_gpp
                    + self.active_power_loss_in_transf_gpp
                )
                * 10 ** -3
                * self.line_length
            )
        )
        u4 = 17 * sqrt(
            (self.line_length / 16)
            + (
                self.calculated_active_power_on_gpp
                + self.active_power_loss_in_transf_gpp
            )
            * 10 ** -3
        )

        sum_u = u1 + u2 + u3 + u4
        self.supply_line_voltage = sum_u / 4

    def calc_of_additional_power(
        self, main_workshop: Workshop, additional_workshops: [Workshop]
    ):
        main_workshop.calculated_reactive_power_сумма = (
            main_workshop.calculated_reactive_power
        )
        main_workshop.calculated_active_power_summ = (
            main_workshop.calculated_active_power
        )

        for дополнительный_цех in additional_workshops:
            main_workshop.calculated_reactive_power_сумма += (
                дополнительный_цех.calculated_reactive_power
            )
            main_workshop.calculated_active_power_summ += (
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
    expenses = 0

    def expenses__calculation(self):
        for i in range(workshops_counts):
            self.к_л += self.workshops[i].line_investments
            self.к_эа += self.workshops[i].capital_investments_switcher
            self.к_тп += self.workshops[i].transformer_price
            self.с_ал += self.workshops[i].deprecation_value
            self.с_пт += self.workshops[i].cost_losses_of_electric_energy_in_tp
            self.с_пл += self.workshops[i].electricity_loss_cost
            self.с_аэа += self.workshops[i].deprecation_value_for_switcher
            self.с_атп += self.workshops[i].deprecation_value_for_tp

    def expenses_calculation(self):
        self.expenses = (
            0.15 * (self.к_л + self.к_эа + self.к_тп)
            + self.с_ал
            + self.с_пт
            + self.с_пл
            + self.с_аэа
            + self.с_атп
        )

    def __init__(self):
        self.workshops = []
        for i in range(len(enterprise_powers)):
            workshop = Workshop(
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
            self.workshops.append(workshop)

        # print(self.workshops[0].calculated_active_power  , 'calculated_active_power')
        self.calc_of_active_power_with_loss_account()
        self.calc_of_reactive_power_with_loss_account()
        # print(self.active_power_with_loss_account, 'active_power_with_loss_account')
        # print(self.reactive_power_with_loss_account, 'reactive_power_with_loss_account')
        self.calc_of_compensating_devices()
        self.calc_of_power_loss_in_compensational_devices()
        self.calc_of_calculated_active_power_on_gpp()
        self.calc_of_calculated_reactive_power_on_tires_gpp()
        self.calc_of_calculated_full_power_on_tires_gpp()
        self.calc_of_loss_active_power_in_transf_gpp()
        self.calc_of_loss_reactive_power_in_transf_gpp()
        self.calc_of_full_calculated_load_with_account_of_power_losses()
        self.calc_of_supply_line_voltage()
        self.expenses__calculation()
        self.expenses_calculation()

        # учет РП(мощности тп + мощности тп)

        for _, (main_workshop_index, additional_workshops_indexes) in enumerate(
            additional_workshops_map.items()
        ):
            main_workshop = self.workshops[main_workshop_index]
            additional_workshops = []

            if additional_workshops_indexes[0] == no_additional_workshops:
                main_workshop.calculated_reactive_power_сумма = (
                    main_workshop.calculated_reactive_power
                )
                main_workshop.calculated_active_power_summ = (
                    main_workshop.calculated_active_power
                )
                continue

            for additional_workshop_index in additional_workshops_indexes:
                additional_workshops.append(self.workshops[additional_workshop_index])

            self.calc_of_additional_power(
                main_workshop=main_workshop, additional_workshops=additional_workshops
            )


# noinspection NonAsciiCharacters
def main():
    workshops_conatiner = WorkshopsContainer()

    rp1 = "1 rp"
    rp2 = "2 rp"
    transf1 = "1 тр"
    transf2 = "2 тр"

    # категории = [2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 3, 3, 2, 2, 2]

    schemas_variants = [
        [rp1, transf1, transf2, rp2] if workshop_category == 3 else [transf2]
        for workshop_category in workshops_categories
    ]

    for index, workshop_category in enumerate(workshops_categories):
        if workshop_category == 3:
            if (
                workshops_conatiner.workshops[index].full_calculated_transformers_load
                <= 330
            ):
                pass

            elif (
                workshops_conatiner.workshops[index].full_calculated_transformers_load
                > 330
                and workshops_conatiner.workshops[
                    index
                ].full_calculated_transformers_load
                <= 660
            ):

                schemas_variants[index].remove(рп1)

            elif (
                workshops_conatiner.workshops[index].full_calculated_transformers_load
                > 660
            ):
                schemas_variants[index].remove(рп1)

                schemas_variants[index].remove(рп2)

    # print(schemas_variants)
    # ------------------------s
    list_of_all_schemas_variants_with = list(
        itertools.product(*schemas_variants)
    )
    workshops_counts = 15
    # ------------------------e

    # ------------------------s
    # arr = []

    # for schema_variant in list_of_all_schemas_variants_with:
    #     # print(schema_variant, 'schema_variant')
    #     zeros_and_ones = [0] * workshops_counts
    #     for workshop in range(len(schema_variant)):
    #         if schema_variant[workshop] == тр2 or schema_variant[workshop] == тр1:
    #             zeros_and_ones[workshop] = 1
    #         else:
    #             zeros_and_ones[workshop] = 0
    #     arr.append(zeros_and_ones)
    # ------------------------e

    # print(arr[2], 'arr[2]')
    # array.append(zeros_and_ones)
    # print(arr[2])

    # ------------------------s
    # for schema_variant in arr:
    #     for workshop in range(len(schema_variant)):
    #         if schema_variant[workshop] == 1:
    #             schema_variant[workshop] = [-1]
    #         elif schema_variant[workshop] == 0:
    #             schema_variant[workshop] = []

    #     for workshop in range(len(schema_variant)):
    #         if schema_variant[workshop] == []:
    #             for workshop_selected in range(len(schema_variant)):
    #                 if schema_variant[workshop_selected] == [-1]:
    #                     schema_variant[workshop].append(workshop_selected)
    # ------------------------e

    # print(arr[200], 'arr 200')
    # print(list_of_all_schemas_variants_with[500], 'индекс 500')
    # print(arr[500], 'индекс 500')
    # print(arr[20], 'индекс 20')
    # print(list_of_all_schemas_variants_with[20], 'индекс 20')
    # print(arr[2], 'index 1')
    # print(arr[5000], 'index 2')
    # print(zeros_and_ones)
    # ar = [[[-1], [-1], [-1], [-1], [-1], [-1], [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 14], [-1], [-1], [-1], [-1], [-1], [-1], [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 14], [-1], [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 14]],
    # [[-1], [-1], [-1], [-1], [-1], [-1], [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11], [-1], [-1], [-1], [-1], [-1], [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11], [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11], [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11], [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11]]]

    # ------------------------s
    # all_possible_connections_list = []
    # for список in arr:
    #     all_possible_connections_list.append(
    #         list(itertools.product(*список)))
    # print(all_possible_connections_list)
    # ------------------------e

    # all_possible_connections_list = list(itertools.product(*ar[0]))
    # print(all_possible_connections_list[1000])
    # print(all_possible_connections_list)
    # result = []
    # for p in all_possible_connections_list:
    #     result.append(p)
    # print(result)
    # print(len(result))
    # # for key, value in карта_доп_цехов:

    # варианты_подключений = []

    # for index, workshop_category in enumerate(workshops_categories):
    #     if workshop_category == 3:
    #         for i in range(len(workshops_categories)):
    #             вариант_подключения = {}
    #             for key, val in вариант_подключения:
    #                 вариант_подключения[key] = [index]

    # варианты_подключений[i] =

    # print(workshops_conatiner.power_loss_in_compensational_devices)
    # print(workshops_conatiner.full_calculated_power_on_tires_gpp)
    # print(workshops_conatiner.full_calculated_load_with_account_of_power_losses)
    # print(workshops_conatiner.supply_line_voltage, 'supply_line_voltage')
    # print(workshops_conatiner.workshops[0].calculated_active_power, 'calculated_active_power workshop 1')
    # print(workshops_conatiner.workshops[0].calculated_reactive_power, 'calculated_reactive_power workshop 1')
    # print(workshops_conatiner.workshops[0].calculated_active_power_summ, 'calculated_active_power_summ workshop 1')
    # print(workshops_conatiner.workshops[0].calculated_reactive_power_сумма, 'calculated_reactive_power_сумма workshop 1')
    # print(workshops_conatiner.workshops[3].calculated_active_power, 'calculated_active_power workshop 4')
    # print(workshops_conatiner.workshops[6].calculated_active_power, 'calculated_active_power workshop 7 ')
    # print(workshops_conatiner.workshops[9].calculated_active_power, 'calculated_active_power workshop 10 ')
    # print(workshops_conatiner.workshops[3].calculated_active_power_summ, 'calculated_active_power_summ')
    # print(workshops_conatiner.workshops[0].transformer_load_coeff, 'transformer_load_coeff')
    # print(workshops_conatiner.workshops[0].reactive_power_loss_ts,
    #       'reactive_power_loss_ts')
    # print(workshops_conatiner.workshops[0].calculated_load_a_gpp_tp,
    #       'calculated_load_a_gpp_tp')
    # print(workshops_conatiner.workshops[0].calculated_load_n_gpp_tp,
    #       'calculated_load_n_gpp_tp')
    # print(workshops_conatiner.workshops[1].calculated_load_a_gpp_tp, 'calculated_load_a_gpp_tp')
    # print(workshops_conatiner.workshops[1].calculated_load_n_gpp_tp, 'calculated_load_n_gpp_tp')
    # print(workshops_conatiner.workshops[0].calculated_current_load_n_gpp_tp, 'calculated_current_load_n_gpp_tp')
    # print(workshops_conatiner.workshops[1].permissible_voltage_loss, 'permissible_voltage_loss')
    # print(workshops_conatiner.workshops[1].permissible_economic_current_density_loss, 'permissible_economic_current_density_loss')
    # print(workshops_conatiner.workshops[0].permissible_current_load_a_gpp_tp, 'permissible_current_load_a_gpp_tp')
    # print(workshops_conatiner.workshops[6].full_calculated_transformers_load, 'full_calculated_transformers_load')
    # print(workshops_conatiner.workshops[3].calculated_active_power_summ, 'calculated_active_power_summ')
    # print(workshops_conatiner.workshops[3].calculated_active_power, 'calculated_active_power')
    # print(workshops_conatiner.workshops[6].calculated_active_power, 'calculated_active_power')
    # print(workshops_conatiner.workshops[0].transformer_load_coeff, 'transformer_load_coeff')
    print(
        workshops_conatiner.workshops[0].electricity_loss_cost,
        "electricity_loss_cost",
    )
    # print(workshops_conatiner.workshops[0].cost_losses_of_electric_energy_in_tp, 'cost_losses_of_electric_energy_in_tp')
    print(
        workshops_conatiner.workshops[11].transformer_load_coeff,
        "transformer_load_coeff",
    )
    print(
        workshops_conatiner.workshops[11].full_calculated_transformers_load,
        "full_calculated_transformers_load",
    )
    print(
        workshops_conatiner.workshops[11].calculated_transformers_power,
        "calculated_transformers_power",
    )
    print(
        workshops_conatiner.workshops[11].preselected_transformer,
        "предварительно_selected_transformer",
    )
    print(workshops_conatiner.workshops[11].selected_transformer, "selected_transformer")
    print(workshops_conatiner.expenses, "expenses")
    # print(workshops_conatiner.workshops[9].length, 'length')
    # print(workshops_conatiner.workshops[1].length, 'length')
    # print(workshops_conatiner.workshops[0].complex_reactive_power_after_compensation, 'complex_reactive_power_after_compensation')

    # print(workshops_conatiner.workshops[1].chosen_cable, 'chosen_cable workshop 2')
    # print(workshops_conatiner.workshops[0].color_material_consumption, 'color_material_consumption')
    # print(workshops_conatiner.workshops[0].cable_current, 'cable_current')

    # print(workshops_conatiner.workshops[0].load_coeff, 'коэффициент_загрузк��')
    # for workshop in workshops_conatiner.workshops:
    #     print(workshop.full_calculated_transformers_load, 'full_calculated_transformers_load ')


if __name__ == "__main__":
    main()
