#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Nomogram-based Total Points Calculation and Disease-Free Survival (DFS) Prediction Using NAR-TRG Scores

This model is derived from a cohort at Sun Yat-sen University Cancer Center and externally validated at Fujian Provincial Cancer Hospital.

Calculation methods for variables:
1) NAR score: Calculated using Valentini's formula based on clinical T stage (cT), pathological T stage after neoadjuvant therapy (ypT), and pathological N stage (pN).
2) TRG score: Modified Mandard grading system ranging from 1-4; grades 1-2 are classified as TRGlow, and grades 3-4 as TRGhigh.
3) CA19-9: Considered positive if greater than 35 U/ml.
4) Tumor differentiation: Categorized as Well, Moderate, or Poor.

The app computes:
1) NAR score and corresponding group (low, medium, high)
2) TRG group (low vs. high)
3) Combined NAR-TRG score (1-3)
4) Nomogram total points and 3-year and 5-year Disease-Free Survival (DFS) predictions.
"""

import math
import streamlit as st

# R-derived parameters
beta_nar = 0.4836
beta_diff = 0.3897
beta_ca19 = 0.5679
LP_base = -0.5138644
S0_3yr = 0.938052
S0_5yr = 0.9052835
LP_max_case = 1.410968
LP_max_rel = LP_max_case - LP_base
POINT_SCALE = 220.0 / LP_max_rel

# Functions for calculating scores and survival
def calculate_nar(cT: int, pT: int, pN: int) -> float:
    numerator = (5*pN - 3*(cT - pT) + 12)
    return (numerator**2) / 9.61

def classify_nar(nar_value: float) -> str:
    if nar_value < 8.0:
        return "NARlow"
    elif nar_value > 16.0:
        return "NARhigh"
    else:
        return "NARmedium"

def classify_trg(trg_value: int) -> str:
    return "TRGlow" if trg_value in [1, 2] else "TRGhigh"

def determine_nar_trg_score(nar_group: str, trg_group: str) -> int:
    if nar_group in ["NARlow", "NARmedium"] and trg_group == "TRGlow":
        return 1
    elif nar_group in ["NARlow", "NARmedium"] and trg_group == "TRGhigh":
        return 2
    else:
        return 3

def predict_points_and_survival(nar_trg: int, poorly_diff: int, ca19_flag: int) -> dict:
    LP_patient = (beta_nar * nar_trg + beta_diff * poorly_diff + beta_ca19 * ca19_flag)
    LP_rel = LP_patient - LP_base
    points = max(LP_rel * POINT_SCALE, 0.0)
    surv_3yr = S0_3yr ** math.exp(LP_rel)
    surv_5yr = S0_5yr ** math.exp(LP_rel)
    return {"Points": points, "DFS_3yr": surv_3yr, "DFS_5yr": surv_5yr}

# Streamlit Interface
def main():
    st.title("Nomogram-based Total Points and DFS Prediction Using NAR-TRG Scores")

    st.markdown("""
    **Instructions:**
    - Select clinical T stage (cT), pathological T stage after neoadjuvant therapy (ypT), and pathological N stage (pN).
    - Choose Modified Mandard TRG grade (1-4).
    - Select tumor differentiation status.
    - Input the numeric value of CA19-9.
    """)

    cT_choice = st.selectbox("Clinical T Stage (cT)", [1,2,3,4], index=0)
    pT_choice = st.selectbox("Pathological T Stage (ypT)", [0,1,2,3,4], index=0)
    pN_choice = st.selectbox("Pathological N Stage (pN)", [0,1,2], index=0)

    nar_score = calculate_nar(cT_choice, pT_choice, pN_choice)
    nar_group = classify_nar(nar_score)
    st.write(f"**NAR Score**: {nar_score:.2f} ({nar_group})")

    trg_choice = st.selectbox("Modified Mandard TRG (1-4)", [1,2,3,4], index=2)
    trg_group = classify_trg(trg_choice)
    st.write(f"**TRG Group**: {trg_group}")

    nar_trg_combined = determine_nar_trg_score(nar_group, trg_group)
    st.write(f"**Combined NAR-TRG Score**: {nar_trg_combined}")

    diff_choice = st.selectbox("Tumor Differentiation", ["Well", "Moderate", "Poor"], index=1)
    diff_is_poor = (diff_choice == "Poor")

    ca19_input = st.number_input("CA19-9 Value (U/ml)", min_value=0.0, value=10.0, step=1.0)
    ca19_over35 = (ca19_input > 35.0)

    if st.button("Compute Total Points and 3-year/5-year DFS"):
        results = predict_points_and_survival(nar_trg_combined, diff_is_poor, ca19_over35)
        st.write(f"**Nomogram Total Points**: {results['Points']:.2f} / 220")
        st.write(f"**Predicted 3-year DFS**: {results['DFS_3yr']*100:.1f}%")
        st.write(f"**Predicted 5-year DFS**: {results['DFS_5yr']*100:.1f}%")

if __name__ == "__main__":
    main()
