#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
A simple interactive web application (using Streamlit) for calculating:

1) The Neoadjuvant Rectal (NAR) score based on cT, ypT, and pN.
2) The NAR-TRG classification (1, 2, or 3) based on NAR (<8, 8-16, >16) and TRG (1-2 vs. 3-4).
3) CA19-9 as a continuous variable input, then converted to a binary variable (>35 or not).
4) Tumor differentiation (well, moderate, poor).
5) A Cox-model-based total points and 3-year/5-year survival probabilities.

How to run:
    1) Install streamlit if not already:  pip install streamlit
    2) Save this script as, e.g.,  "nar_nomogram_app.py"
    3) In a terminal, run:  streamlit run nar_nomogram_app.py
    4) A browser tab will open, allowing interactive input.
    
Note:
    - The numeric values for the Cox model coefficients, baseline survival, etc.,
      come from prior analyses in R (rms package). These are for DEMO.
    - Adjust them to match your actual fitted model.
"""

import math
import streamlit as st

# -------------------------------------------------------------------
# 1. Define functions for NAR score and classification
# -------------------------------------------------------------------

def calculate_nar(cT: int, pT: int, pN: int) -> float:
    """
    Calculate the Neoadjuvant Rectal (NAR) score using Valentini's formula:
    
    NAR = [5*pN - 3*(cT - pT) + 12]^2 / 9.61
    
    All values (cT, pT, pN) are integers with ranges:
      cT in [1..4], pT in [0..4], pN in [0..2].
    """
    # Ensure nonnegative results in bracket by formula definition
    numerator = (5*pN - 3*(cT - pT) + 12)
    return (numerator**2) / 9.61

def classify_nar(nar_value: float) -> str:
    """
    Classify NAR score into 'low', 'medium', or 'high'.
      - NAR < 8  => low
      - 8..16    => medium
      - >16      => high
    """
    if nar_value < 8.0:
        return "NARlow"
    elif nar_value > 16.0:
        return "NARhigh"
    else:
        return "NARmedium"

def classify_trg(trg_value: int) -> str:
    """
    Mandard TRG classification:
      - TRG 1-2 => TRGlow
      - TRG 3-4 => TRGhigh
    (Some systems use TRG up to 5; typically TRG>=3 is high. 
     Adjust as needed.)
    """
    if 1 <= trg_value <= 2:
        return "TRGlow"
    else:
        return "TRGhigh"

def determine_nar_trg_score(nar_group: str, trg_group: str) -> int:
    """
    Combine NAR group and TRG group into NAR-TRG scores:
        Score 1 => (NARlow TRGlow) or (NARmedium TRGlow)
        Score 2 => (NARlow TRGhigh) or (NARmedium TRGhigh)
        Score 3 => (NARhigh TRGlow) or (NARhigh TRGhigh)
    """
    if nar_group in ["NARlow", "NARmedium"] and trg_group == "TRGlow":
        return 1
    elif nar_group in ["NARlow", "NARmedium"] and trg_group == "TRGhigh":
        return 2
    else:
        # nar_group == "NARhigh"
        return 3

# -------------------------------------------------------------------
# 2. Cox-model logic for total points & 3y/5y survival
#    (Example: using previously derived parameters from an R 'cph' model)
# -------------------------------------------------------------------

# Example log(HR) from cph (DEMO values):
BETA_NARTRG  = {1:0.0, 2:0.5, 3:1.0}  
# ^ Suppose we treat NAR-TRG score (1,2,3) as if: 
#    Score1 => 0.0, Score2 => 0.5, Score3 => 1.0 in log(HR). 
# Just for demonstration. Replace with real estimates 
# or a more complex approach if your final model is different.

BETA_DIFF    = 0.4   # e.g. for 'poor' differentiation vs. else
BETA_CA19    = 0.56  # e.g. for CA19>35 vs. not

# Suppose baseline has (NAR-TRG=1 => 0.0 log(HR)), non-poor => 0, CA19<=35 => 0
LP_BASE      = -0.51 # if the real baseline from R is e.g. -0.51
S0_3YR       = 0.94  # 3-year baseline
S0_5YR       = 0.90  # 5-year baseline

# For points mapping, let's say the highest relative LP => 220 points 
# (Here we do a simplified version. If you want, compute your own 'max scenario'.)
LP_MAX_SCENARIO = 1.41   # e.g. from R
POINT_SCALE     = 220.0 / (LP_MAX_SCENARIO - LP_BASE)

def compute_cox_points_survival(
    nar_trg_score: int,
    diff_is_poor: bool,
    ca19_over_35: bool
) -> (float, float, float):
    """
    Returns (total_points, survival_3yr, survival_5yr).
    Using a simplified Cox approach:
    
      LP(patient) = BETA_NARTRG[nar_trg_score] 
                  + (BETA_DIFF if diff_is_poor else 0)
                  + (BETA_CA19 if ca19_over_35 else 0)
      LP_rel = LP(patient) - LP_BASE
      Points = LP_rel * POINT_SCALE
      Surv(t) = (S0_t) ^ exp(LP_rel)
    """
    lp_patient = BETA_NARTRG.get(nar_trg_score, 0.0)
    if diff_is_poor:
        lp_patient += BETA_DIFF
    if ca19_over_35:
        lp_patient += BETA_CA19

    lp_rel = lp_patient - LP_BASE
    points = lp_rel * POINT_SCALE
    if points < 0:
        points = 0.0

    # 3y / 5y survival
    surv_3y = S0_3YR ** (math.exp(lp_rel))
    surv_5y = S0_5YR ** (math.exp(lp_rel))
    return (points, surv_3y, surv_5y)

# -------------------------------------------------------------------
# 3. Build a Streamlit UI
# -------------------------------------------------------------------
def main():
    st.title("NAR Score & Nomogram Calculator (Demo)")

    st.markdown("""
    **Instructions**:  
    - Select cT, ypT, pN from drop-downs.  
    - Choose Mandard TRG.  
    - Input tumor differentiation (well, moderate, or poor).  
    - Input numeric CA19-9, we'll check if >35.  
    - We compute:  
      1) NAR score, NAR group (low, medium, high)  
      2) TRG group (low vs. high)  
      3) NAR-TRG score (1,2,3)  
      4) Nomogram total points & 3y/5y survival (demo).
    """)
    
    # cT, ypT, pN (all are ints)
    cT_choice  = st.selectbox("cT (1-4)", options=[1,2,3,4], index=0)
    pT_choice  = st.selectbox("ypT (0-4)", options=[0,1,2,3,4], index=0)
    pN_choice  = st.selectbox("pN (0-2)", options=[0,1,2], index=0)

    # Compute NAR
    nar_score = calculate_nar(cT_choice, pT_choice, pN_choice)
    nar_group = classify_nar(nar_score)

    st.write(f"**NAR score** = {nar_score:.2f}, classified as **{nar_group}**")

    # TRG (Mandard)
    trg_choice = st.selectbox("TRG (1-5)", options=[1,2,3,4,5], index=2)
    trg_group  = classify_trg(trg_choice)
    st.write(f"TRG = {trg_choice}, => {trg_group}")

    # NAR-TRG
    nar_trg_combined = determine_nar_trg_score(nar_group, trg_group)
    st.write(f"NAR-TRG score = **{nar_trg_combined}**")

    # Tumor differentiation
    diff_choice = st.selectbox("Tumor differentiation", 
        options=["Well", "Moderate", "Poor"], index=1)

    diff_is_poor = (diff_choice == "Poor")

    # CA19-9 input
    ca19_input = st.number_input("CA19-9 Value (U/ml)", min_value=0.0, value=10.0, step=1.0)
    ca19_over35 = (ca19_input > 35.0)

    # Press button to compute final
    if st.button("Compute Nomogram Points & Survival"):
        # In this DEMO, we treat the NAR-TRG score as an ordinal factor => 1,2,3
        # Then apply our simplified model:
        points, surv3, surv5 = compute_cox_points_survival(
            nar_trg_combined,
            diff_is_poor,
            ca19_over35
        )
        
        st.write(f"**Nomogram total points** = {points:.2f} (max ~220)")
        st.write(f"**3-year survival** = {surv3*100:.1f}%")
        st.write(f"**5-year survival** = {surv5*100:.1f}%")

if __name__ == "__main__":
    main()
