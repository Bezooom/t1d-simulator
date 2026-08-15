# Release notes — `v1.1.0`

**Tag:** `v1.1.0`  
**Date:** 2026-08-15  
**License:** MIT  
**Repo:** https://github.com/Bezooom/t1d-simulator  

---

## English

### t1d_simulator v1.1.0 — Verified implant geometry

What changed since v1.0.0: the main claim is no longer “we have a simulator”, but a **checked rule for implant shape and size**.

**Not a clinical product.** In silico oxygen bounds + literature check.

#### Geometry (the science)

- Sufficient Krogh length \(L^*\), compared to the MM BVP (mean \(L_\mathrm{num}/L^* = 1.08\), **0** certificate violations)
- **Forbidden:** planar SQ half-thickness ≥ 250 µm at ≥ 50 M/ml; 400 µm capsules at 50 M/ml (matches Papas / Papabathini viability drops)
- **Allowed:** omental sphere ~200 µm / 80 M/ml; thin SQ planar ≤ 100 µm at ≤ 50 M/ml
- Edmonton-like dose (70 kg, 11k IEQ/kg) → **~650 cm²** thin omental sheet, not a thicker slab
- Native 150 µm IEQ already ~38% anoxic core in SQ; pseudoislets **≤ 120 µm**
- Avascular wait: thick SQ keeps ~14% mass after 18 human days; certified omental sphere \(f_\mathrm{anox}=0\)

Reproduce:

```bash
python3 t1d_simulator/design_bounds.py
python3 t1d_simulator/dose_transport.py
python3 t1d_simulator/scaleup_science.py
```

#### Product

- `screen_design --sweep` + heatmaps (`reports/sweeps/`)
- Preprint **v1.3** EN + RU
- GitHub language switch: [English](https://github.com/Bezooom/t1d-simulator/blob/master/README.md) · [Русский](https://github.com/Bezooom/t1d-simulator/blob/master/README.ru.md)
- Tests: **44 / 44** (`CUDA_VISIBLE_DEVICES="" python3 t1d_simulator/verify_model.py`)

#### Links

| | |
|--|--|
| Russian landing | https://github.com/Bezooom/t1d-simulator/blob/master/README.ru.md |
| Preprint EN | https://github.com/Bezooom/t1d-simulator/blob/master/docs/manuscript_biorxiv_en.md |
| Preprint RU | https://github.com/Bezooom/t1d-simulator/blob/master/docs/manuscript_biorxiv_ru.md |
| Clinician one-pager (RU) | https://github.com/Bezooom/t1d-simulator/blob/master/docs/endocrinologist_brief_ru.md |
| \(L^*\) report | https://github.com/Bezooom/t1d-simulator/blob/master/reports/design_bounds/verified_design_space.md |

```
git clone https://github.com/Bezooom/t1d-simulator.git
cd t1d-simulator
pip install -r t1d_simulator/requirements.txt
```

---

## Русский

### t1d_simulator v1.1.0 — обоснованная геометрия имплантата

С v1.0.0 главное не «есть симулятор», а **проверенное правило формы и размера**.

Это не лекарство и не протокол пересадки.

#### Геометрия

- Достаточная граница Крога \(L^*\), сверка с численным BVP (среднее отношение 1.08, **0** нарушений сертификата)
- **Нельзя:** толстый SQ ≥ 250 мкм при ≥ 50 млн/мл; капсулы 400 мкм при 50 млн/мл (как у Papas / Papabathini)
- **Можно:** сфера в сальнике ~200 мкм / 80 млн/мл; тонкий SQ ≤ 100 мкм при ≤ 50 млн/мл
- Доза как в Эдмонтоне → **~650 см²** тонкого сальникового листа (с A4), не толще блин
- Островок 150 мкм под кожей уже с ~38% мёртвого ядра; псевдоостровки **≤ 120 мкм**
- Пока нет сосудов: толстый SQ к 18 суткам человека ~14% массы; сертифицированная сфера в сальнике не платит этот налог

#### Продукт

- `screen_design --sweep` и тепловые карты
- Препринт **v1.3** на русском и английском
- На GitHub сверху: **English** / **Русский**
- Тесты **44 / 44**

Коротко для коллег-эндокринологов: [пояснение](https://github.com/Bezooom/t1d-simulator/blob/master/docs/endocrinologist_brief_ru.md)
