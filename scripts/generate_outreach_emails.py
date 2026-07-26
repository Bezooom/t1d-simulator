# -*- coding: utf-8 -*-
"""
Automated Outreach Email Generator
Generates personalized, non-spammy academic PI and biotech R&D outreach drafts.
"""

import os
import sys

TARGET_CONTACTS = [
    {
        "name": "Dr. David Mooney",
        "institution": "Harvard SEAS",
        "type": "Academic PI",
        "focus": "biomaterial hydrogel macropores and FBR reduction",
        "offer": "We would be glad to run a free 1D/3D multiphysics oxygen transport simulation tailored to your specific alginate pore dimensions and cell loading density."
    },
    {
        "name": "Dr. Daniel Anderson",
        "institution": "MIT Koch Institute",
        "type": "Academic PI",
        "focus": "zwitterionic hydrogel modifications for immune evasion",
        "offer": "We can evaluate your novel polymer SMILES structures using our GNN anti-fibrotic model to predict fibrotic thickness (L_fib) and oxygen flux in vivo."
    },
    {
        "name": "R&D Team",
        "institution": "Sana Biotechnology",
        "type": "Biotech",
        "focus": "hypoimmune iPSC-derived cell therapies",
        "offer": "Our digital twin engine models population survival under B2M-/- / CD47+ edit sets, 0-48h IBMIR kinetics, and omental pouch retrievability."
    },
    {
        "name": "Cell Therapy Team",
        "institution": "Seraxis / Sernova",
        "type": "Biotech",
        "focus": "pre-vascularized pouches and micro-organoid constructs",
        "offer": "We offer a quantitative comparison of graft oxygenation, VEGF angiogenesis feedback, and steatosis risks across omental vs. subcutaneous sites."
    }
]

EMAIL_TEMPLATE = """Subject: In Silico Failure Screening (`screen_design`) & Collaboration for {institution}

Dear {name},

I have been following your impressive research on {focus}.

We have developed an open-source, literature-calibrated digital twin engine (`t1d_simulator`) designed to evaluate multiphysics failure modes in stem cell-derived beta-cell therapies—specifically oxygen diffusion limits (Krogh limit), FBR fibrotic capsule growth, and 0–48h acute IBMIR kinetics.

Using our new automated construct screening tool (`screen_design`), {offer}

Our preprint manuscript, open-science package, and screening tools are accessible here:
- Executable Screening Tool: `python3 t1d_simulator/screen_design.py`
- Manuscript: docs/manuscript_biorxiv_en.md (bioRxiv PDF / Zenodo DOI metadata)
- Decision Matrix: docs/site_comparison_matrix.md

We would welcome the opportunity to run a zero-cost computational failure screening for your group's custom geometry or discuss model co-development.

Best regards,

Pavel V. Naumov
Bioengineering Research Group | t1d_simulator
naumov122@gmail.com
"""

def generate_emails():
    output_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "reports", "outreach_emails")
    os.makedirs(output_dir, exist_ok=True)
    
    print("=== Generating Personalized Outreach Email Drafts ===")
    generated_files = []
    
    for idx, c in enumerate(TARGET_CONTACTS, start=1):
        content = EMAIL_TEMPLATE.format(
            name=c["name"],
            institution=c["institution"],
            focus=c["focus"],
            offer=c["offer"]
        )
        safe_inst = c['institution'].replace(' ', '_').replace('/', '_').lower()
        filename = f"email_{idx}_{safe_inst}.txt"
        filepath = os.path.join(output_dir, filename)
        
        with open(filepath, "w", encoding="utf-8") as f:
            f.write(content)
            
        generated_files.append(filepath)
        print(f"  [{idx}] Generated draft for {c['name']} ({c['institution']}) -> {filename}")
        
    print(f"\n[SUCCESS] {len(generated_files)} email drafts created in reports/outreach_emails/")
    return generated_files

if __name__ == "__main__":
    generate_emails()
