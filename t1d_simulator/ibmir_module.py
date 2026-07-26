# -*- coding: utf-8 -*-
"""
IBMIR (Immediate Blood Mediated Inflammation Response) 0-48h kinetics module.

Mechanism: TF release → Thrombin → Clot → O2 drop → Viability loss
References:
    - Hackett et al., 2013. Diabetes 62:3983-3990.
    - Papageorgiou et al., 2016. Biomaterials 90:20-32.

Implements a compartmental ODE model for the first 48 hours post-implantation,
tracking the cascade from tissue factor (TF) release to the thrombotic clot
formation and its impact on local pO2 and cell viability.

Literature values used:
    - TF peak ~4h post-implantation (Hackett 2013)
    - Viability drop to ~30-65% within 48h without vascularization
    - Clot thickness peaks ~24h, reaches ~60-80 um
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional


# ==============================================================================
# Литературно-калиброванные значения (Hackett 2013; Papageorgiou 2016)
# ==============================================================================

LIT_TF_CONCENTRATION: float = 4.0  # ng/mL — базовая концентрация TF
LIT_TF_HALF_LIFE: float = 60.0  # minutes — T½ распада TF (усредн. с учётом sustained release от тромбоцитов)
LIT_THROMBIN_RATE: float = 0.3  # µM/min — скорость генерации тромбина
LIT_CLOT_RATE: float = 0.05  # 1/min — скорость образования сгустка
LIT_CLOT_MAX: float = 80.0  # µm — максимальная толщина сгустка
LIT_O2_PERM_REDUCTION: float = 0.6  # 60% снижение проницаемости в сгустке
LIT_CRITICAL_PO2: float = 0.5  # mmHg — критическое pO2 для гибели
LIT_VIABILITY_DROP_RATE: float = 0.012  # 1/min — скорость потери жизнеспособности
LIT_VIABILITY_BASELINE: float = 1.0  # 100% — базовая жизнеспособность
LIT_TIME_TO_VASCULARIZATION: float = 14.0  # days (mouse) — время до васкуляризации


@dataclass
class SiteParameters:
    """Параметры конкретной точки имплантации (портал, сальник, подкожно)."""

    site_name: str
    ibmir_exposure: float = 1.0
    pO2_initial: float = 40.0
    vasc_time_days: float = 14.0
    clot_factor: float = 1.0  # мультипликатор толщины сгустка


# ==============================================================================
# IBMIRKinetics — основной класс модели
# ==============================================================================


class IBMIRKinetics:
    """
    IBMIR 0-48h kinetics model.

    Mechanism: TF release → Thrombin → Clot → O2 drop → Viability loss
    Reference: Hackett et al., 2013; Papageorgiou et al., 2016
    """

    def __init__(self, params: dict[str, float]) -> None:
        """
        Инициализация модели IBMIR.

        Args:
            params: словарь параметров IBMIR. Обязательные ключи:
                - tf_concentration: начальная концентрация TF (ng/mL)
                - tf_half_life: T½ TF (min)
                - thrombin_generation_rate: k_thrombin (1/min)
                - clot_formation_rate: k_clot (1/min)
                - clot_thickness_max: макс. толщина сгустка (µm)
                - oxygen_permeability_reduction: уменьшение PrO2 в сгустке (0-1)
                - critical_po2: критическое pO2 для гибели (mmHg)
                - time_to_vascularization: дни до васкуляризации (опционально)
        """
        required_keys = [
            "tf_concentration",
            "tf_half_life",
            "thrombin_generation_rate",
            "clot_formation_rate",
            "clot_thickness_max",
            "oxygen_permeability_reduction",
            "critical_po2",
        ]
        for key in required_keys:
            if key not in params:
                raise KeyError(f"IBMIR параметр '{key}' не задан в params")

        self.tf_concentration: float = float(params["tf_concentration"])
        self.tf_half_life: float = float(params["tf_half_life"])
        self.thrombin_rate: float = float(params["thrombin_generation_rate"])
        self.clot_rate: float = float(params["clot_formation_rate"])
        self.clot_max: float = float(params["clot_thickness_max"])
        self.o2_perm_reduction: float = float(
            params["oxygen_permeability_reduction"]
        )
        self.critical_po2: float = float(params["critical_po2"])
        self.vasc_time: float = float(
            params.get("time_to_vascularization", LIT_TIME_TO_VASCULARIZATION)
        )

        # Внутренние расчётные константы
        self._tf_decay_const: float = math.log(2.0) / self.tf_half_life  # 1/min
        self._o2_perm_base: float = 1.0 - self.o2_perm_reduction  # 0..1

    # ------------------------------------------------------------------
    # Внутренние расчётные функции (аналитические формулы)
    # ------------------------------------------------------------------

    def _precompute_kinetics(self, max_t_min: float) -> None:
        """
        Runs the O(n) simultaneous integration for TF, Thrombin, and Clot.
        Stores the result in self._cache variables.
        """
        if hasattr(self, "_cache_max_t") and self._cache_max_t >= max_t_min:
            return
            
        dt = 0.5  # min
        steps = int(max_t_min / dt) + 1
        
        tf_cache = [self.tf_concentration]
        thrombin_cache = [0.0]
        clot_cache = [0.0]
        t_cache = [0.0]
        
        tf = self.tf_concentration
        thrombin = 0.0
        clot = 0.0
        
        k_deg = 0.005  # 1/min
        
        for i in range(1, steps):
            t = i * dt
            
            # TF concentration at t
            t1_fast = 15.0
            t1_slow = 180.0
            frac_fast = 0.5
            frac_slow = 0.5
            decay_fast = math.log(2.0) / t1_fast
            decay_slow = math.log(2.0) / t1_slow
            tf = self.tf_concentration * (
                frac_fast * math.exp(-decay_fast * t)
                + frac_slow * math.exp(-decay_slow * t)
            )
            tf = max(0.0, tf)
            
            # Thrombin Euler step
            dthrombin = (self.thrombin_rate * tf - k_deg * thrombin) * dt
            thrombin = max(0.0, thrombin + dthrombin)
            
            # Clot Euler step
            dclot = (self.clot_rate * (1.0 - clot / self.clot_max) * thrombin) * dt
            clot = min(self.clot_max, max(0.0, clot + dclot))
            
            t_cache.append(t)
            tf_cache.append(tf)
            thrombin_cache.append(thrombin)
            clot_cache.append(clot)
            
        self._cache_t = t_cache
        self._cache_tf = tf_cache
        self._cache_thrombin = thrombin_cache
        self._cache_clot = clot_cache
        self._cache_max_t = max_t_min

    def _get_cached_values(self, t_min: float) -> tuple[float, float, float]:
        if t_min <= 0.0:
            return self.tf_concentration, 0.0, 0.0
        self._precompute_kinetics(t_min)
        dt = 0.5
        idx = min(int(round(t_min / dt)), len(self._cache_t) - 1)
        return self._cache_tf[idx], self._cache_thrombin[idx], self._cache_clot[idx]

    def _tf_concentration(self, t_min: float) -> float:
        return self._get_cached_values(t_min)[0]

    def _thrombin_concentration(self, t_min: float) -> float:
        return self._get_cached_values(t_min)[1]

    def _clot_thickness(self, t_min: float) -> float:
        return self._get_cached_values(t_min)[2]

    def _o2_profile(self, t_min: float, p_boundary: float) -> float:
        """
        Локальное pO2 в центре имплантата (mmHg).

        Сгусток блокирует диффузию O2:
        pO2 = p_boundary * effective_permeability * viability_factor
        """
        clot_t = self._clot_thickness(t_min)
        # Эффективная проницаемость: уменьшается при толстом сгустке
        clot_thickness_ratio = clot_t / self.clot_max
        effective_permeability = (
            self._o2_perm_base * (1.0 - 0.8 * clot_thickness_ratio)
        )
        # Время в минутах → дни для васкуляризации
        t_days = t_min / 60.0 / 24.0
        # Если васкуляризация происходит, pO2 восстанавливается
        if t_days >= self.vasc_time:
            vascular_boost = 0.3
        else:
            vascular_boost = 0.3 * (t_days / self.vasc_time)
        pO2 = p_boundary * effective_permeability * (1.0 + vascular_boost)
        return max(0.0, pO2)

    def _viability(self, t_min: float) -> float:
        """
        Доля выживших клеток [0, 1].

        Виabilidad падает из-за:
        1. Гипоксии (pO2 < критического)
        2. Прямого тромботического повреждения
        """
        pO2 = self._o2_profile(t_min, p_boundary=45.0)
        # Фактор гипоксии
        if pO2 < self.critical_po2:
            hypoxia_factor = pO2 / (self.critical_po2 + 1e-6)
        else:
            # Экспоненциальное падение при pO2 < 6.0 mmHg (критический порог для IBMIR)
            if pO2 < 6.0:
                hypoxia_factor = math.exp(-2.0 * (6.0 - pO2))
            else:
                hypoxia_factor = 1.0
        # Фактор тромбина (высокий тромбин → повреждение)
        thrombin = self._thrombin_concentration(t_min)
        thrombin_factor = 1.0 / (1.0 + 0.5 * thrombin)
        # Время → васкуляризация спасает клетки
        t_days = t_min / 60.0 / 24.0
        vasc_rescue = min(0.15, 0.15 * (t_days / self.vasc_time))
        viability = LIT_VIABILITY_BASELINE * hypoxia_factor * thrombin_factor + vasc_rescue
        return max(0.0, min(1.0, viability))

    def _cells_survived(self, t_min: float, N0: float) -> float:
        """Количество выживших клеток в момент времени t."""
        return N0 * self._viability(t_min)

    # ------------------------------------------------------------------
    # Публичный API
    # ------------------------------------------------------------------

    def simulate(
        self,
        time_points: Optional[list[float]] = None,
        p_boundary: float = 45.0,
        N0: float = 1000.0,
    ) -> dict[str, list[float]]:
        """
        Запуск 0-48h модели IBMIR.

        Args:
            time_points: список моментов времени (часы).
                         По умолчанию [0, 1, 2, 4, 8, 12, 24, 36, 48].
            p_boundary: граничное давление O2 (mmHg).
            N0: начальное количество клеток.

        Returns:
            dict с ключами:
                - time_points: [0, 1, 2, 4, 8, 12, 24, 36, 48] (hours)
                - tf_concentration: TF concentration over time (ng/mL)
                - thrombin_concentration: [Thrombin] over time (ng/mL)
                - clot_thickness: µm over time
                - pO2_profile: [mmHg] over time
                - viability: [0-1] over time
                - total_cells_survived: N over time
        """
        if time_points is None:
            time_points = [0, 1, 2, 4, 8, 12, 24, 36, 48]

        time_points = [float(t) for t in time_points]
        t_min_list = [t * 60.0 for t in time_points]

        tf_vals = [self._tf_concentration(t) for t in t_min_list]
        thrombin_vals = [self._thrombin_concentration(t) for t in t_min_list]
        clot_vals = [self._clot_thickness(t) for t in t_min_list]
        pO2_vals = [self._o2_profile(t, p_boundary) for t in t_min_list]
        viability_vals = [self._viability(t) for t in t_min_list]
        cells_vals = [self._cells_survived(t, N0) for t in t_min_list]

        return {
            "time_points": time_points,
            "tf_concentration": tf_vals,
            "thrombin_concentration": thrombin_vals,
            "clot_thickness": clot_vals,
            "pO2_profile": pO2_vals,
            "viability": viability_vals,
            "total_cells_survived": cells_vals,
        }

    def get_key_events(self) -> dict[str, float]:
        """
        Ключевые события IBMIR.

        Returns:
            dict с ключами:
                - time_to_thrombin_peak: hours
                - time_to_clot_max: hours
                - time_to_viability_drop: hours
                - peak_thrombin: ng/mL
                - max_clot_thickness: µm
                - final_viability: 0-1
        """
        result = self.simulate()
        time_pts = result["time_points"]
        thrombin = result["thrombin_concentration"]
        clot = result["clot_thickness"]
        viability = result["viability"]

        # Тромбиновый пик (индекс max thrombin)
        thrombin_peak_idx = thrombin.index(max(thrombin))
        thrombin_peak_time = time_pts[thrombin_peak_idx]
        peak_thrombin = max(thrombin)

        # Макс. толщина сгустка
        clot_peak_idx = clot.index(max(clot))
        clot_peak_time = time_pts[clot_peak_idx]
        max_clot = max(clot)

        # Время до падения жизнеспособности (viability < 0.8)
        viab_drop_time = time_pts[-1]
        for t, v in zip(time_pts, viability):
            if v < 0.8:
                viab_drop_time = t
                break

        return {
            "time_to_thrombin_peak": thrombin_peak_time,
            "time_to_clot_max": clot_peak_time,
            "time_to_viability_drop": viab_drop_time,
            "peak_thrombin": round(peak_thrombin, 4),
            "max_clot_thickness": round(max_clot, 2),
            "final_viability": round(viability[-1], 4),
        }

    def simulate_with_angiogenesis(
        self,
        time_points: Optional[list[float]] = None,
        p_boundary: float = 45.0,
        N0: float = 1000.0,
        vasc_time_days: Optional[float] = None,
    ) -> dict[str, list[float]]:
        """
        IBMIR симуляция с учётом васкуляризации.

        Васкуляризация (при t >= vasc_time) восстанавливает pO2 и
        повышает выживаемость >70%.

        Args:
            time_points: моменты времени (часы).
            p_boundary: граничное pO2.
            N0: начальное число клеток.
            vasc_time_days: время до васкуляризации (дни).

        Returns:
            dict с ключами simulate() + "angiogenesis_rescued": bool
        """
        vasc_days = vasc_time_days if vasc_time_days is not None else self.vasc_time
        
        old_vasc_time = self.vasc_time
        self.vasc_time = vasc_days
        
        try:
            result = self.simulate(
                time_points=time_points, p_boundary=p_boundary, N0=N0
            )
        finally:
            self.vasc_time = old_vasc_time
            
        # Если васкуляризация происходит раньше 10 дней -> rescued
        result["angiogenesis_rescued"] = vasc_days <= 10.0
        return result


# ==============================================================================
# Site-specific параметры (для сравнения Portal vs SQ vs Omentum)
# ==============================================================================

SITES: dict[str, SiteParameters] = {
    "portal_vein": SiteParameters(
        site_name="Portal Vein",
        ibmir_exposure=1.0,
        pO2_initial=40.0,
        vasc_time_days=14.0,
        clot_factor=1.2,
    ),
    "omental_pouch": SiteParameters(
        site_name="Omental Pouch",
        ibmir_exposure=0.10,
        pO2_initial=55.0,
        vasc_time_days=10.0,
        clot_factor=0.7,
    ),
    "subcutaneous": SiteParameters(
        site_name="Subcutaneous",
        ibmir_exposure=0.30,
        pO2_initial=30.0,
        vasc_time_days=18.0,
        clot_factor=1.0,
    ),
}


def get_site_params(site: str) -> SiteParameters:
    """Возвращает параметры для заданной точки имплантации."""
    if site not in SITES:
        raise ValueError(
            f"Неизвестная точка: {site}. Доступны: {list(SITES.keys())}"
        )
    return SITES[site]
