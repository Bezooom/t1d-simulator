# Дорожная карта: цифровой двойник терапии СД1

**Проект:** In silico β-cell / hypoimmune organoid digital twin (`t1d_simulator`)  
**Версия:** 3.0  
**Дата:** 2026-07-26  
**Статус:** living document — **M0–M3 (in-repo) закрыты 100%**; внешние действия + M4–M5 впереди  

**Связанные документы**

| Документ | Роль |
|----------|------|
| `Development_Roadmap.md` | Архив НИОКР-фаз 1–5 (симулятор) |
| `Development_Plan_Research_Impact.md` | Стратегия impact (архив идей) |
| `Experimental_Biochemist_Plan.md` | In vitro / in vivo — только с партнёром |
| `docs/wet_lab_validation_protocol.md` | SOP валидации endpoints |
| `docs/manuscript_biorxiv_en.md` | Preprint EN (submit-ready) |
| `docs/grant_proposal_breakthrough_t1d.md` | Черновик grant aims (~$1.25M) |
| `Critical_Analysis_*.md` | Limitations → methods preprint |
| `zenodo.json` | Метаданные citeable software release |

---

## 0. Executive status (2026-07-26)

### 0.1. Вехи in-repo

| Веха | Описание | Статус |
|------|----------|--------|
| **M0 Engineering** | `parameters.yaml`, `literature_params.yaml`, `param_loader.py`, EN README, CI | ✅ **100%** |
| **M1 Calibration** | O₂ PDE / Krogh; Papas RMSE **10.42%**; Papabathini RMSE **11.42%**; site matrix; IBMIR 0–48 h | ✅ **100%** |
| **M2 Open Science** | Preprint EN/RU + PDF; tornado `uncertainty_analysis.py`; CLI `screen_design.py`; `zenodo.json` | ✅ **100%** |
| **M3 Outreach Ready** | Web screening в UI; `run_partner_screening.py`; hybrid emails в `reports/outreach_emails/` | ✅ **100%** |
| **Тест-сюита** | `verify_model.py` | ✅ **40 / 40 green** |
| **M4 Joint science** | Co-authorship / grant / model v2 по partner data | ⬜ внешний + joint |
| **M5 Wet-lab** | Партнёр + Experimental Plan E1+ | ⬜ только с lab |

> **Уточнение терминов:** M3 здесь = **Outreach Ready** (пакет, скрипты, письма, UI готовы).  
> **M3 Signal** (внешний KPI: ≥2 meaningful scientific replies) наступает **после** External Actions ниже.

### 0.2. Что остаётся вне репозитория (External Actions)

| # | Действие | Критерий done | Блокирует |
|---|----------|---------------|-----------|
| **E1** | Релизный тег **`v1.0.0`** на GitHub + публикация кода | Tag + public repo; Zenodo DOI через `zenodo.json` / GitHub–Zenodo | Citeability, письма |
| **E2** | Подача preprint на **bioRxiv** (Synthetic Biology / Bioengineering) | Preprint ID / URL | Outreach credibility |
| **E3** | Рассылка черновиков из `reports/outreach_emails/` | 4+ отправленных; log дат | M3 Signal → M4 |
| **E4** | (После ответов) grant Breakthrough T1D / joint figure | Submission / MoU | M4–M5 |

**Порядок:** E1 → E2 → E3 (письма с **живыми** URL demo/DOI/preprint).  
M4–M5 **не** стартуют «автоматически из кода» — они стартуют после **ответов партнёров** и/или подачи grant.

---

## 1. Северная звезда

### 1.1. Миссия (12–24 месяца)

Сделать **калиброванный, воспроизводимый, открытый digital twin** дизайнов замещения β-клеток: failure modes макроинкапсуляции + сценарии hypoimmune organoids — так, чтобы lab/biotech могли **скринить геометрию, плотность, сайт и edit-set до wet-lab** и получать **честные negative predictions**.

### 1.2. Научное обещание

> Макроинкапсуляция имеет предсказуемые multiphysics failure modes.  
> Digital twin отсекает невозможные density / geometry / site и ранжирует hypoimmune + omentum сценарии **до** животного.

**Формула позиционирования:**  
*independent multiphysics failure analysis + design-space screening; offering model co-development.*

### 1.3. Принцип рычага

| Уровень | Влияние | Статус |
|---------|---------|--------|
| **A.** Open tools + models | Ускоряет чужие дизайны | ✅ in-repo |
| **B.** Validation + citeable release | Доверие | 🟡 артефакты ✅; **DOI/preprint URL = E1–E2** |
| **C.** Co-authorship / partner | Протокол wet-lab | ⬜ M4 |
| **D.** Свой клеточный продукт | Годы / капитал | ⬜ не цель соло |

---

## 2. Закрытые deliverables M0–M3 (реестр)

### 2.1. M0 — Engineering ✅

| Артефакт | Путь |
|----------|------|
| Runtime parameters | `t1d_simulator/parameters.yaml` |
| Literature-calibrated params | `data/literature_params.yaml` |
| Loader | `t1d_simulator/param_loader.py`, `config/loader.py` |
| EN README / packaging | `t1d_simulator/README.md`, `pyproject.toml` |
| CI skeleton | `.github/workflows/ci.yml` |

### 2.2. M1 — Calibration ✅

| Артефакт | Путь / метрика |
|----------|----------------|
| O₂ PDE + multiphysics | `simulator.py` |
| IBMIR 0–48 h | `ibmir_module.py`, `organoid_simulator.py` |
| Site comparison | `docs/site_comparison_matrix.md` + tests |
| Krogh / VEGF / IBMIR benchmarks | `reports/benchmarks/` |
| Papas et al. viability | **RMSE 10.42%** |
| Papabathini et al. viability | **RMSE 11.42%** |
| One-command reproduce | `reports/benchmarks/reproduce_benchmarks.py` |

### 2.3. M2 — Open Science (in-repo) ✅

| Артефакт | Путь |
|----------|------|
| Preprint EN | `docs/manuscript_biorxiv_en.md` **v1.1** (submit-ready draft) |
| Preprint RU | `docs/manuscript_biorxiv_ru.md` **v1.1** (PDF/HTML могут отставать — пересобрать перед подачей) |
| Figures / graphical abstract | `docs/figures/` |
| Tornado / uncertainty | `t1d_simulator/uncertainty_analysis.py` |
| Design-space CLI | `t1d_simulator/screen_design.py` |
| Zenodo metadata | `zenodo.json` (`v1.0.0`) |
| Export demo pack | `t1d_simulator/export_demo_pack.py` |

> **M2 Public** (bioRxiv live + Zenodo DOI) = External Actions **E1–E2**.

### 2.4. M3 — Outreach Ready ✅

| Артефакт | Путь |
|----------|------|
| Interactive UI screening | `t1d_simulator/app.py`, `ui/panels/` |
| Partner batch screening | `scripts/run_partner_screening.py` |
| 1-pager / brief / contacts | `docs/outreach_*.md` |
| Hybrid email drafts | `reports/outreach_emails/email_1…4_*.txt` |
| MoU template | `docs/collab_agreement_mou_template.md` |
| Email generator | `scripts/generate_outreach_emails.py` |

**Стартовые адресаты (черновики готовы):**

1. Harvard SEAS — Mooney Lab (`email_1_harvard_seas.txt`)  
2. MIT Koch — Anderson Lab (`email_2_mit_koch_institute.txt`)  
3. Sana Biotechnology (`email_3_sana_biotechnology.txt`)  
4. Seraxis / Sernova (`email_4_seraxis___sernova.txt`)  

> **M3 Signal** = отправка (E3) + ≥2 meaningful replies → вход в M4.

### 2.5. Quality gate

| Проверка | Результат |
|----------|-----------|
| `python3 t1d_simulator/verify_model.py` | **40 / 40 green** |
| Core physics + PINN + GNN + IBMIR + sites | covered |
| `screen_design` + tornado sensitivity | tests 39–40 |
| Open-science pack integrity | validate in suite |

---

## 3. Архитектура дальше

```
[M0–M3 IN-REPO ✅] ──► E1 GitHub v1.0.0 + Zenodo DOI
                              │
                              ▼
                       E2 bioRxiv submit
                              │
                              ▼
                       E3 Outreach send ──► M3 Signal (≥2 replies)
                              │
              ┌───────────────┴───────────────┐
              ▼                               ▼
        M4 Joint science              Grant Breakthrough T1D
     (joint figure / MoU / v2)         (blueprint ~$1.25M)
              │                               │
              └──────────► M5 Wet-lab E1+ ◄───┘
                     (только с партнёром)
```

---

## 4. External Actions — runbook

### E1. GitHub release `v1.0.0` + Zenodo DOI

**Цель:** citeable software artifact.

1. Убедиться, что remote public (или станет public при release).  
2. В корне: **`LICENSE`** (MIT) — ✅ `LICENSE`; текст релиза — ✅ `docs/RELEASE_NOTES_v1.0.0.md`.  
3. Не коммитить: `venv/`, `.env`, secrets, крупные кэши.  
4. Tag:

```bash
git add -A
git status   # review
git commit -m "Release v1.0.0: M0–M3 digital twin pack"
git tag -a v1.0.0 -m "t1d_simulator v1.0.0 — multiphysics T1D digital twin"
git push origin master
git push origin v1.0.0
```

5. GitHub → **Create release** from tag `v1.0.0` — вставить body из `docs/RELEASE_NOTES_v1.0.0.md` (блок под `---`).  
6. Zenodo GitHub integration **или** manual upload с `zenodo.json` → записать **DOI** сюда:

| Поле | Значение |
|------|----------|
| GitHub release | `⬜` URL |
| Zenodo DOI | `⬜` 10.5281/zenodo.XXXX |

### E2. bioRxiv preprint

**Цель:** публичный scientific record.

1. PDF: **`docs/manuscript_biorxiv_en.pdf`** (primary, v1.1) или `docs/manuscript_biorxiv_ru.pdf`. Пересборка: `python3 scripts/build_formatted_manuscript.py en ru`.  
2. Category: **Synthetic Biology** / **Bioengineering**.  
3. Code availability: GitHub URL + Zenodo DOI (после E1).  
4. Data/params: `data/literature_params.yaml`, `reports/benchmarks/`.  
5. Записать:

| Поле | Значение |
|------|----------|
| bioRxiv ID / URL | `⬜` |
| Posted date | `⬜` |

### E3. Outreach send (M3 → Signal)

**Цель:** scientific conversations, не vanity opens.

1. Подставить в письма: **GitHub**, **demo URL** (если deployed), **Zenodo DOI**, **bioRxiv URL**.  
2. Отправить минимум 4 черновика; расширить до ~15 PI (см. `docs/outreach_contacts_list.md`).  
3. Offer: free `run_partner_screening` / `screen_design` на *их* geometry.  
4. Log:

| Дата | Адресат | Канал | Ответ | Next |
|------|---------|-------|-------|------|
| | Mooney Lab | email | ⬜ | |
| | Anderson Lab | email | ⬜ | |
| | Sana | email | ⬜ | |
| | Seraxis/Sernova | email | ⬜ | |

**M3 Signal done when:** ≥2 meaningful replies **или** joint figure discussion.

### E4. Grant (параллельно после E2 + ≥1 signal)

- Blueprint: `docs/grant_proposal_breakthrough_t1d.md` (~$1.25M, 3 years).  
- Не подавать «в вакуум» без preprint URL и хотя бы soft interest lab — слабее review.  
- Target: Breakthrough T1D Strategic Research / uni IT+biophysics / NIDDK-style aims.

---

## 5. M4 — Joint science (после M3 Signal)

| # | Задача | Выход | P |
|---|--------|-------|---|
| M4.1 | Joint figure (predicted vs their / public wet data) | shared figure + methods | P0 |
| M4.2 | Model v2 calibration on partner endpoints | release notes + params changelog | P0 |
| M4.3 | MoU / collab agreement | signed from template | P1 |
| M4.4 | Second preprint **или** grant submission | ID / receipt | P0 |
| M4.5 | Optional: public Streamlit/HF demo always-on | URL | P1 |

**Exit M4:** co-authorship path **или** grant submitted/scored **или** signed MoU + shared data plan.

### Предлагаемый joint figure #1 (минимальная стоимость)

Predicted vs measured **core viability** at 2–3 thicknesses/densities under controlled hypoxia (MIN6 / pseudoislets, 48 h) — без CRISPR и без large animal.

---

## 6. M5 — Wet-lab (только с партнёром)

Источник: `Experimental_Biochemist_Plan.md`, `docs/wet_lab_validation_protocol.md`.

| Фаза | Содержание | In silico support |
|------|------------|-------------------|
| **E1** | 3D print gyroid + hypoxia assay | STL + predicted viable fraction |
| **E2** | SBAA/CBAA + FBR histology | GNN ranking + L_fib |
| **E3** | VEGF + STZ mouse | VEGF + death window |
| **E4** | anti-IBMIR / omental pouch / iCasp9 checklist | organoid modules |

**Gate:** письменный интерес + protocol owner + ethics/IACUC на стороне lab.

---

## 7. Near-term трек (опционально, не бренд)

| Направление | Когда |
|-------------|-------|
| AID (`aid_controller.py`) / CGM hypo-risk | После E2, spare bandwidth |
| Track1/2 repurposing → short note | P2 |

Не смешивать с hypoimmune organoid pitch в одном cold email.

---

## 8. Метрики

### 8.1. Уже достигнуто (in-repo)

| Метрика | Факт |
|---------|------|
| Literature params | ✅ YAML + sources |
| Hard benchmarks | ✅ Papas 10.42%, Papabathini 11.42% |
| Preprint draft EN/RU/PDF | ✅ |
| Screen + uncertainty modules | ✅ |
| Outreach pack | ✅ 4 emails + 1-pager |
| Tests | ✅ 40/40 |

### 8.2. Следующие 90 дней (external + M4)

| Метрика | Target |
|---------|--------|
| GitHub `v1.0.0` + Zenodo DOI | 1 |
| bioRxiv live | 1 |
| Emails sent | ≥10–15 (min 4) |
| Meaningful replies | ≥2 |
| Partner free screens delivered | ≥1–2 |
| Grant submitted | 0–1 (stretch) |

---

## 9. Анти-план

1. Новые «фазы биомиметики» вместо E1–E3.  
2. Overclaim: TIR/CGM как clinical trial; GNN N≈46 как validated product.  
3. Cold email **без** preprint/DOI URL.  
4. Соло wet-lab / animals без партнёра.  
5. Ждать «идеальный 3D CFD» перед bioRxiv.  
6. Vanity metrics (opens) вместо scientific replies.

---

## 10. Definition of Done

**Code (уже норма для M0–M3):**

- [x] Module + assert in `verify_model`  
- [x] Params from YAML  
- [x] Limitations reflected in Critical Analysis / manuscript  

**External:**

- [ ] E1 DOI recorded in this ROADMAP §4  
- [ ] E2 preprint URL recorded  
- [ ] E3 log ≥4 sent; M3 Signal when ≥2 replies  

**M4 task:**

- [ ] Shared artifact (figure / data / MoU) with named external party  

---

## 11. Одностраничный summary

```
M0 Engineering .............. ✅ 100%
M1 Calibration .............. ✅ 100%  (Papas 10.42% / Papabathini 11.42%)
M2 Open Science (in-repo) ... ✅ 100%  (preprint pack, screen, tornado, zenodo.json)
M3 Outreach Ready ........... ✅ 100%  (UI, partner script, 4 emails)
Tests ....................... ✅ 40/40

NEXT (external):
  E1  git tag v1.0.0 → Zenodo DOI
  E2  bioRxiv submit (EN PDF / md package)
  E3  send outreach emails with live URLs
  →  M3 Signal (≥2 replies) → M4 Joint → M5 Wet-lab (partner)
  →  Grant Breakthrough T1D (blueprint ready)
```

**Северная звезда:** ускорять правильные wet-lab решения и резать тупики инкапсуляции — не симулировать «уже cure».

**Главный next step:** **E1 → E2 → E3** (не новый код).

---

## 12. Changelog

| Дата | Версия | Изменение |
|------|--------|-----------|
| 2026-07-23 | 1.0 | Единая карта M0–M5, kickstart |
| 2026-07-26 | 2.0 | Baseline M0 ✅ / M1 ~80%; фокус M2 ship |
| 2026-07-26 | **3.0** | **M0–M3 in-repo 100%**; 40 tests; реестр deliverables; External Actions E1–E4 runbook; M4–M5 как post-signal; M3 = Ready vs Signal |
)
