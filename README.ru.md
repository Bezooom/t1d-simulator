<a id="readme-top"></a>
<!-- best-readme-template -->
<div align="center">

<h1 align="center">Цифровой двойник СД1 и гипоиммунных органоидов</h1>

<p align="center">
  Открытый откалиброванный по литературе мультифизический цифровой двойник для инкапсуляции островков и дизайна гипоиммунных органоидов при СД1.
  <br />
  <br />
  <a href="https://github.com/Bezooom/t1d-simulator/tree/HEAD/docs"><strong>Документация »</strong></a> · <a href="https://github.com/Bezooom/t1d-simulator/issues/new?labels=bug">Сообщить об ошибке</a> · <a href="https://github.com/Bezooom/t1d-simulator/issues/new?labels=enhancement">Предложить улучшение</a> · <a href="README.md">English</a>
</p>

<p align="center">
  <a href="https://github.com/Bezooom/t1d-simulator/actions/workflows/ci.yml"><img src="https://img.shields.io/github/actions/workflow/status/Bezooom/t1d-simulator/ci.yml?style=for-the-badge&label=CI&logo=githubactions&logoColor=white" alt="CI" /></a>
  <a href="https://github.com/Bezooom/t1d-simulator/releases/latest"><img src="https://img.shields.io/github/v/release/Bezooom/t1d-simulator?style=for-the-badge" alt="Release" /></a>
  <a href="https://github.com/Bezooom/t1d-simulator/blob/HEAD/LICENSE"><img src="https://img.shields.io/badge/License-MIT-blue?style=for-the-badge" alt="License" /></a>
  <a href="https://github.com/Bezooom/t1d-simulator/stargazers"><img src="https://img.shields.io/github/stars/Bezooom/t1d-simulator?style=for-the-badge" alt="Stars" /></a>
  <a href="https://github.com/Bezooom/t1d-simulator/commits"><img src="https://img.shields.io/github/last-commit/Bezooom/t1d-simulator?style=for-the-badge" alt="Last commit" /></a>
  <a href="https://github.com/Bezooom/t1d-simulator/issues"><img src="https://img.shields.io/github/issues/Bezooom/t1d-simulator?style=for-the-badge" alt="Issues" /></a>
</p>

</div>

<details open>
<summary><b>📑 Содержание</b></summary>

- [📖 О проекте](#-о-проекте)
- [Для эндокринолога (одна страница)](#для-эндокринолога-одна-страница)
- [Что умеет код](#что-умеет-код)
- [🧱 Технологии (Built With)](#-технологии-built-with)
- [Быстрый старт](#быстрый-старт)
- [Документация](#документация)
- [Проверка на чужих данных](#проверка-на-чужих-данных)
- [Цитирование](#цитирование)
- [🗺️ Дорожная карта](#-дорожная-карта)
- [🤝 Участие в разработке](#-участие-в-разработке)
- [Лицензия](#лицензия)
- [📬 Контакты](#-контакты)
- [🙏 Благодарности](#-благодарности)

</details>

## 📖 О проекте

**Язык / Language:** [English](README.md) · **Русский**

[![Tests](https://github.com/Bezooom/t1d-simulator/actions/workflows/ci.yml/badge.svg)](https://github.com/Bezooom/t1d-simulator/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Препринт RU](https://img.shields.io/badge/препринт-RU-blue.svg)](docs/manuscript_biorxiv_ru.md)
[![Preprint EN](https://img.shields.io/badge/preprint-EN-red.svg)](docs/manuscript_biorxiv_en.md)
[![Release](https://img.shields.io/badge/release-v1.1.0-brightgreen.svg)](https://github.com/Bezooom/t1d-simulator/releases/tag/v1.1.0)

Открытый расчёт **Павла Наумова**: какой формы и толщины можно сажать инкапсулированные $\beta$-клетки при СД1, а какой — нельзя (кислород не доходит до середины). Это не лекарство и не протокол пересадки.

**Главный результат:** граница Крога $L^*$. Толстый подкожный лист ($\ge 250$ мкм полутолщины при $\ge 50$ млн/мл) и капсулы радиусом $400$ мкм на типичной плотности — **выше** $L^*$. Тонкий лист у сальника или мелкие сферы — нет. Доза уровня Эдмонтона — это **площадь** (~650 см², с лист A4), не «ещё потолще».

![Графический абстракт](docs/figures/graphical_abstract.png)

---

## Для эндокринолога (одна страница)

[Пояснение простыми словами](docs/endocrinologist_brief_ru.md) · [статья RU](docs/manuscript_biorxiv_ru.md) · [статья EN](docs/manuscript_biorxiv_en.md)

---

## Что умеет код

- **Кислород (PDE)** — пластина / цилиндр / сфера, кинетика Михаэлиса–Ментен, фиброз, АФК.
- **Граница $L^*$** — аналитика + сверка с численным решателем (`design_bounds.py`).
- **Доза → площадь** — сколько см² сальника на 70 кг (`dose_transport.py`).
- **Островок и окно до сосудов** — ядро нативного IEQ, мышь 8 сут. vs человек 18 сут. (`scaleup_science.py`).
- **Скрининг** — `screen_design` (один дизайн или `--sweep` + heatmap).
- **IBMIR 0–48 ч** — TF → тромбин → сгусток, гепарин.
- **PINN / GNN / UI** — ускорение, разведочный скрининг полимеров, Streamlit. Не главный claim.

---

## 🧱 Технологии (Built With)

<p>
<img src="https://img.shields.io/badge/Python-3-3776AB?style=for-the-badge&logo=python&logoColor=white" alt="Python" />
<img src="https://img.shields.io/badge/NumPy--013243?style=for-the-badge&logo=numpy&logoColor=white" alt="NumPy" />
<img src="https://img.shields.io/badge/SciPy--8CAAE6?style=for-the-badge&logo=scipy&logoColor=white" alt="SciPy" />
<img src="https://img.shields.io/badge/PyTorch--EE4C2C?style=for-the-badge&logo=pytorch&logoColor=white" alt="PyTorch" />
<img src="https://img.shields.io/badge/Streamlit--FF4B4B?style=for-the-badge&logo=streamlit&logoColor=white" alt="Streamlit" />
</p>

## Быстрый старт

```bash
git clone https://github.com/Bezooom/t1d-simulator.git
cd t1d-simulator
pip install -r t1d_simulator/requirements.txt

export CUDA_VISIBLE_DEVICES=""
python3 t1d_simulator/verify_model.py          # 44 теста
streamlit run t1d_simulator/app.py

# один дизайн
python3 t1d_simulator/screen_design.py --geometry spherical --radius 200 --density 80 --site omental_pouch

# сетка + картинка
python3 t1d_simulator/screen_design.py --sweep \
  --geometry planar --site subcutaneous --no-heparin \
  --densities 20,50,80,100 --radii 100,150,200,250 \
  --heatmap reports/sweeps/heatmap_viability_planar_sq.png

# проверенные границы
python3 t1d_simulator/design_bounds.py
python3 t1d_simulator/dose_transport.py
python3 t1d_simulator/scaleup_science.py
```

---

## Документация

| Документ | Русский | English |
|----------|---------|---------|
| Этот README | [README.ru.md](README.ru.md) | [README.md](README.md) |
| Препринт v1.3 | [manuscript_biorxiv_ru.md](docs/manuscript_biorxiv_ru.md) | [manuscript_biorxiv_en.md](docs/manuscript_biorxiv_en.md) |
| Для клинициста | [endocrinologist_brief_ru.md](docs/endocrinologist_brief_ru.md) | — |
| Граница $L^*$ | [verified_design_space.md](reports/design_bounds/verified_design_space.md) | то же |
| Доза → площадь | [dose_vs_krogh.md](reports/dose_transport/dose_vs_krogh.md) | то же |
| Островок + сосуды | [remaining_biotech_bounds.md](reports/scaleup/remaining_biotech_bounds.md) | то же |
| Кому писать | [outreach_playbook.md](docs/outreach_playbook.md) | — |
| Оглавление docs | [docs/README.md](docs/README.md) | [docs/README.md](docs/README.md) |

---

## Проверка на чужих данных

| Бенчмарк | Статья | Итог |
|----------|--------|------|
| Предел Крога | Secomb et al., 2004 | PASS |
| Сферы альгината | Papas et al., 2007 | RMSE **10.42%** |
| Планарный макродевайс | Papabathini et al., 2023 | RMSE **11.42%** |

---

## Цитирование

```bibtex
@software{t1d_digital_twin_2026,
  title = {T1D Digital Twin \& Hypoimmune Organoid Simulator},
  author = {Naumov, Pavel V.},
  year = {2026},
  url = {https://github.com/Bezooom/t1d-simulator},
  version = {1.0.0},
  license = {MIT}
}
```

## 🗺️ Дорожная карта

Подробный план: [ROADMAP.md](ROADMAP.md). Полный список задач и идей — в [issues](https://github.com/Bezooom/t1d-simulator/issues).

## 🤝 Участие в разработке

Вклад в проект приветствуется. Если есть идея, как сделать лучше, откройте issue или pull request:

1. Сделайте форк репозитория
2. Создайте ветку (`git checkout -b feature/my-idea`)
3. Сделайте коммит (`git commit -m 'feat: my idea'`)
4. Отправьте ветку (`git push origin feature/my-idea`)
5. Откройте pull request

## Лицензия

MIT — файл [LICENSE](LICENSE). Вопросы: `naumov122@gmail.com`.

## 📬 Контакты

**Павел Наумов** (Bezooom) — [naumovai.ru](https://naumovai.ru) · [X @Bezooooom](https://x.com/Bezooooom) · naumov122@gmail.com

Проект: [https://github.com/Bezooom/t1d-simulator](https://github.com/Bezooom/t1d-simulator)

## 🙏 Благодарности

- Опубликованные физиологические параметры, на которые ссылаются [документация](docs) и препринт
- Метаданные [Zenodo](https://zenodo.org) (`zenodo.json`)
- [Best-README-Template](https://github.com/othneildrew/Best-README-Template)
- [Shields.io](https://shields.io)

<p align="right">(<a href="#readme-top">наверх</a>)</p>
