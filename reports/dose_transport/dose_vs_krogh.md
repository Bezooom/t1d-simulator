# Dose vs oxygen: why encapsulation trials run out of area

**Date:** 2026-08-15  
**Patient:** 70 kg · Edmonton-like **11000 IEQ/kg**  
**Mass:** 770000 IEQ = **1201 million cells** (1560 cells/IEQ)  
**Oxygen:** human $V_{\max}=1.5\times 10^{-16}$ mol/(cell·s), 2% alginate, certificate $0.9\,L^*$  
**Two $p_{\mathrm{crit}}$:** survival $0.5$ mmHg vs **function / GSIS $5$ mmHg** ($K_{M,\mathrm{ins}}$)

## The biotech conflict

To replace insulin you need a **fixed cell mass**.  
Krogh forbids packing that mass into a thick lump.  
So the only remaining degree of freedom is **surgical area** (omental sheet, pouch, or many microcapsules).

$$
V = \frac{N_{\mathrm{cells}}}{\rho},\qquad
A_{\mathrm{planar,\ 2-sided}} = \frac{V}{2\,L^{\ast}_{\mathrm{func}}}
$$

$L^{\ast}_{\mathrm{func}}$ is the same zero-order bound as before, but with $p_c=5$ mmHg. It is **strictly smaller** than the survival bound ($p_c=0.5$). A graft can be “alive” on paper and still not secrete.

## Numbers that matter

| Scenario | $L^*_{\mathrm{surv}}$ | $L^*_{\mathrm{func}}$ | Graft volume | Surgical implication |
|----------|---------------------------|----------------------------|--------------|----------------------|
| Omentum planar 80 M/ml | 135.1 µm | 129.4 µm | 15.02 ml | **644.5 cm²** two-sided sheet |
| SQ planar 50 M/ml | 125.7 µm | 115.8 µm | 24.02 ml | **1153.0 cm²** two-sided (worse if one-sided: 2306.0 cm²) |
| Omentum sphere 80 M/ml | 234.0 µm | 224.2 µm | 15.02 ml | **436.0** spheres at $R=0.9 L^*_{\mathrm{func}}$ |

A credit-card is ~46 cm². A **~644 cm²** omental sheet is large but surgically thinkable.  
A **~2306 cm²** one-sided SQ plaque is why thick subcutaneous macrodevices fail clinically: people thicken the slab instead of growing area, and $L$ crosses $L^*$.

## Literature-typical devices vs $L^*$

| device | L_or_R_um | L_star_surv_um | L_star_func_um | oxygen_survival | oxygen_function | ibmir_ok | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Alginate microcapsule (Papas-like) | 400.0 | 217.8 | 200.5 | False | False | True | FAIL |
| Planar SQ macrodevice, 500 µm full (2-sided → L=250) | 250.0 | 125.7 | 115.8 | False | False | True | FAIL |
| Planar SQ thin, 200 µm full (L=100) | 100.0 | 125.7 | 115.8 | True | True | True | PASS |
| Omental sphere (this study, recommended) | 200.0 | 234.0 | 224.2 | True | True | True | PASS |
| Portal naked / uncoated sphere | 150.0 | 199.2 | 187.5 | True | True | False | FAIL |

Reading: Papas-like **400 µm** capsules at 50 M/ml SQ sit **above** functional $L^*$ (~190 µm class) — core hypoxia is not a surprise, it is the bound. The recommended omental 200 µm / 80 M/ml sphere is the first row that is oxygen-functional *and* not portal-naked.

## Mouse vs human (why papers look better than trials)

| | Mouse | Human |
|--|-------|-------|
| $V_{\max}$ used here | $1.2\times 10^{-16}$ | $1.5\times 10^{-16}$ (we certified on human) |
| Time to vessels | ~7–10 d | ~14–21 d |
| Device scale | mm | must cover Edmonton mass → tens of cm² |

$L^*\propto 1/\sqrt{Q}$. Human OCR + longer avascular window + 1000× more cells is a **quantitative** reason encapsulation “works in mice” and dies in SQ human devices. Not a new slogan — the same formula.

## FBR (order-of-magnitude)

A 50 µm fibrotic shell does not change $Q$, it drops boundary $pO_2$. Cutting $p_b$ from 30 → 20 mmHg shrinks planar $L^*_{\mathrm{func}}$ by $\sqrt{(20-5)/(30-5)}=\sqrt{0.6}\approx 0.77$. Certified thickness falls ~23%; required area rises ~30%. FBR is an **area tax**, not a separate mystery.

## What a lab can use tomorrow

1. Do not design planar SQ thicker than **~$L^*_{\mathrm{func}}$** at the intended density (50 M/ml SQ ≈ **115.8 µm** half-thickness before safety factor).  
2. If you need Edmonton mass, budget **area first** (table), then print/pouch that area.  
3. Portal oxygen can pass $L^*$ and still lose the graft to IBMIR — heparin or leave the blood.  
4. Report **functional** core $pO_2$ (GSIS), not only live/dead.

## Reproduce

```bash
python3 t1d_simulator/dose_transport.py
```

Oxygen certificate details: `reports/design_bounds/verified_design_space.md`.
