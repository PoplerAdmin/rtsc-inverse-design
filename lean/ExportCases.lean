/-
ExportCases: evaluates `TargetLogic.decideT` (the SAME generic definition the theorems are about,
instantiated at `K = ℚ`, exact arithmetic) on a fixed grid of inputs and writes
`export/decide_cases.json` (relative to the current directory) for testing the Python mirror.

Run (from the build directory):  lake env lean ExportCases.lean
-/
import TargetLogic

open TargetLogic

namespace ExportCases

/-! ## JSON helpers (rationals are exported as exact strings "p/q" or "p") -/

def jq (x : ℚ) : String := "\"" ++ toString x ++ "\""
def js (s : String) : String := "\"" ++ s ++ "\""
def jb (b : Bool) : String := if b then "true" else "false"
def jn (n : ℕ) : String := toString n
def jarr (l : List String) : String := "[" ++ ", ".intercalate l ++ "]"
def jobj (l : List (String × String)) : String :=
  "{" ++ ", ".intercalate (l.map fun (k, v) => js k ++ ": " ++ v) ++ "}"
def jopt {α : Type} (f : α → String) : Option α → String
  | some a => f a
  | none => "null"

def tierS : Tier → String
  | .fullDFPT => "fullDFPT" | .coarseDFPT => "coarseDFPT" | .ml => "ML" | .external => "external"
def jprov (p : Prov) : String := jobj [("tier", js (tierS p.tier)), ("run", jn p.run)]

def ruleS : Rule → String
  | .kappa => "N1-KAPPA" | .rowmax => "N1-ROWMAX" | .dos => "N1-DOS" | .iso => "N1-ISO"
  | .nmCert => "NM-CERT" | .nmG => "NM-G"

def assS : Assumption → String
  | .meModel => "ME-MODEL" | .certH => "CERT-H" | .envSound => "ENV-SOUND"
  | .mDiscrete => "M-DISCRETE" | .mSheetIso => "M-SHEET-ISO" | .sheetsCorrect => "SHEETS-CORRECT"
  | .mIso => "M-ISO" | .certNM => "CERT-NM"

def inS : InputName → String
  | .isoA => "A" | .isoLam => "lambda" | .isoM => "m" | .eps => "eps" | .mat => "calA"
  | .dos => "w"

def outLean : Outcome → String
  | .rejectedProven _ => "rejectedProven" | .screenedOut => "screenedOut"
  | .undecided => "undecided" | .pass => "pass"
def outDoc : Outcome → String
  | .rejectedProven _ => "REJECTED_PROVEN" | .screenedOut => "UNDECIDED"
  | .undecided => "UNDECIDED" | .pass => "PASS"
def outRule : Outcome → Option Rule
  | .rejectedProven r => some r
  | _ => none

def jtag : InputTag ℚ → String
  | .fileIncomplete => jobj [("kind", js "fileIncomplete")]
  | .fileCompleteUnverified => jobj [("kind", js "fileCompleteUnverified")]
  | .dropNonpos m => jobj [("kind", js "dropNonpos"), ("massK", jq m)]
  | .negClip m => jobj [("kind", js "negClip"), ("massK", jq m)]
  | .sigmaExcluded r c => jobj [("kind", js "sigmaExcluded"), ("blocksRejection", jb r),
      ("blocksCertification", jb c)]
  | .informational s => jobj [("kind", js "informational"), ("name", js s)]
  | .unknown s => jobj [("kind", js "unknown"), ("name", js s)]

def jvec {n : ℕ} (f : Fin n → ℚ) : String := jarr ((List.ofFn f).map jq)

def jmat (md : MatData ℚ) : String :=
  jobj [("B", jn md.B), ("calA", jarr (List.ofFn fun i => jvec (md.calA i))),
    ("n", jvec md.n), ("v", jvec md.v), ("d_calA", jq md.dCalA), ("prov", jprov md.prov)]

def jdos (dd : DosData ℚ) : String :=
  jobj [("n", jn dd.n), ("w", jvec dd.w), ("d_w", jq dd.dW), ("prov", jprov dd.prov)]

def jinputs (x : Inputs ℚ) : String :=
  jobj [
    ("validity", jobj [("E_F", jopt jq x.val.EF), ("omega_max", jq x.val.omegaMax),
      ("imag_abs_max_K", jq x.val.imagAbsMaxK), ("contains_H_Li_Be", jb x.val.containsHLiBe),
      ("anharmonic_check", jb x.val.anharmonicCheck),
      ("anharm_shift_abs", jq x.val.anharmShiftAbs),
      ("alpha2F_recomputed_anharmonic", jb x.val.alpha2FRecomputedAnharmonic),
      ("mu_star_negative_suspected", jb x.val.muStarNegativeSuspected),
      ("admitted_smearings", jn x.val.admittedSmearings), ("q_meshes", jn x.val.qMeshes),
      ("n_sheets", jn x.val.nSheets)]),
    ("iso", jobj [("A_nom", jq x.iso.A), ("lambda_nom", jq x.iso.lam), ("m_nom", jq x.iso.m),
      ("d_A", jq x.iso.dA), ("d_lambda", jq x.iso.dLam), ("d_m", jq x.iso.dM),
      ("prov_A", jprov x.iso.provA), ("prov_lambda", jprov x.iso.provLam),
      ("prov_m", jprov x.iso.provM)]),
    ("eps", jobj [("eps_sys", jq x.eps.epsSys), ("eps_omega", jq x.eps.epsOmega),
      ("qspread_omega_log", jq x.eps.qSpread), ("prov", jprov x.eps.prov)]),
    ("mat", jopt jmat x.mat), ("dos", jopt jdos x.dos), ("tags", jarr (x.tags.map jtag))]

def jcfg (c : Config ℚ) : String :=
  jobj [("a0", jq c.a0), ("Tc_target", jq c.Tstar),
    ("g_table", jarr (c.tbl.map fun p => jarr [jq p.1, jq p.2])), ("g_inf", jq c.gInf)]

def jcase (cfg : Config ℚ) (id : String) (x : Inputs ℚ) : String :=
  let o := decideT cfg x
  jobj [
    ("id", js id),
    ("inputs", jinputs x),
    ("expected", jobj [
      ("lean_outcome", js (outLean o)), ("status", js (outDoc o)),
      ("rule", jopt (fun r => js (ruleS r)) (outRule o)),
      ("assumptions", jarr (o.assumptions.map (js ∘ assS))),
      ("inputs_used", jarr (o.inputsUsed.map fun n =>
        jobj [("name", js (inS n)), ("prov", jopt jprov (provOf x n))])),
      ("would_fire_ranking_only", jb (o == .screenedOut))]),
    ("gates", jobj [("well_formed", jb (wf x)), ("g_valid", jb (gvalid x)),
      ("blocked", jb (blocked x))]),
    ("rules", jobj (allRules.map fun r =>
      (ruleS r, jobj [("floor", jb (floorR cfg x r)), ("prov_ok", jb (provOK x r))]))),
    ("derived", jobj [("A_hi", jq (Ahi x)), ("lambda_hi", jq (lamHi x)), ("m_hi", jq (mHi x)),
      ("floor_A", jq cfg.floor), ("nm_cert_floor", jq (4 * cfg.floor ^ 2 / lamHi x)),
      ("g_up", jq (cfg.gUp (lamHi x))),
      ("kappa_hi", jopt (fun md => jq (kappaHi md x.eps)) x.mat),
      ("rowmax_hi", jopt (fun md => jq (rowMaxHi md x.eps)) x.mat),
      ("wmin_lo", jopt (fun dd => jq (wminLo dd)) x.dos),
      ("moments_consistent_nominal", jb (decide (x.iso.A ^ 2 ≤ x.iso.lam * x.iso.m / 4)))])]

/-! ## Configuration: doc constants and the certified §0.3 table (`NMTable`, from the certificate JSONs) -/

def cfg : Config ℚ := Config.ofCert (433148145 / 100000000) 300

/-! ## Grid -/

def pf (t : Tier) (r : ℕ) : Prov := ⟨t, r⟩
def full1 : Prov := pf .fullDFPT 1

def val0 : Validity ℚ :=
  { EF := some 100000, omegaMax := 1000, imagAbsMaxK := 0, containsHLiBe := false,
    anharmonicCheck := false, anharmShiftAbs := 0, alpha2FRecomputedAnharmonic := false,
    muStarNegativeSuspected := false, admittedSmearings := 3, qMeshes := 2, nSheets := 1 }

def isoD (A lam m : ℚ) (tA : Tier := .fullDFPT) : IsoData ℚ :=
  ⟨A, lam, m, 1/20, 1/20, 1/20, pf tA 1, full1, full1⟩

def eps0 : EpsData ℚ := ⟨1/10, 1/50, 0, full1⟩

def matSmall (p : Prov) : MatData ℚ :=
  ⟨2, ![![600, 100], ![100, 600]], ![1, 1], ![1, 1], 1/20, p⟩
def matLarge (p : Prov) : MatData ℚ :=
  ⟨2, ![![1500, 100], ![100, 1500]], ![1, 1], ![1, 1], 1/20, p⟩
def dos2 (p : Prov) : DosData ℚ := ⟨2, ![1/2, 1/2], 1/20, p⟩
def dos1 (p : Prov) : DosData ℚ := ⟨1, ![1], 1/20, p⟩

-- Rejected control (A_nom = 500, λ = 1, m = 1e6: consistent, A² = λm/4): N1-ISO fires.
def base : Inputs ℚ := ⟨val0, isoD 500 1 1000000, eps0, none, none, []⟩
/-- Pass control (A_nom = 1250, λ = 1, m = 8e6). -/
def pass0 : Inputs ℚ := ⟨val0, isoD 1250 1 8000000, eps0, none, none, []⟩

def dosBad : DosData ℚ := ⟨2, ![2/5, 1/2], 1/20, full1⟩
def xWsum : Inputs ℚ := ⟨{ val0 with nSheets := 2 }, isoD 1250 1 3000000, eps0, none, some dosBad, []⟩
def xKres : Inputs ℚ := ⟨val0, isoD 1050 1 9000000, eps0, some (matLarge full1), none, []⟩
def xWlow : Inputs ℚ := ⟨{ val0 with nSheets := 2 }, isoD 1250 1 8000000, eps0, none, some ⟨2, ![1/2, 101/200], 0, full1⟩, []⟩
def xLamCoarse : Inputs ℚ := ⟨val0, { (isoD 1150 3 2000000) with provLam := pf .coarseDFPT 1, provM := pf .external 1 }, eps0, none, none, []⟩

def matOpts : List (String × Option (MatData ℚ)) :=
  [("none", none), ("small_full_r1", some (matSmall full1)),
   ("large_full_r1", some (matLarge full1)), ("small_ML", some (matSmall (pf .ml 1))),
   ("small_full_r2", some (matSmall (pf .fullDFPT 2)))]

def gridCases : List (String × Inputs ℚ) := Id.run do
  let mut out := []
  for A in [(500 : ℚ), 1000, 1250, 2000] do
    for m in [(1000000 : ℚ), 3000000, 8000000] do
      for nS in [1, 2] do
        for (mn, mo) in matOpts do
          for dn in [false, true] do
            for tA in [Tier.fullDFPT, Tier.coarseDFPT] do
              let dosO := if dn then some (if nS = 1 then dos1 full1 else dos2 full1) else none
              let lam : ℚ := (([1, 5/2, 8, 20] : List ℚ).find? (fun l => decide (A ^ 2 ≤ l * m / 4))).getD 20
              let x : Inputs ℚ := ⟨{ val0 with nSheets := nS }, isoD A lam m tA, eps0, mo, dosO, []⟩
              let id := s!"grid/A={A}/m={m}/lambda={lam}/nS={nS}/mat={mn}/dos={dn}/tierA={tierS tA}"
              out := out ++ [(id, x)]
  return out

def lamCases : List (String × Inputs ℚ) := Id.run do
  let mut out := []
  for lam in [(1/2 : ℚ), 1, 2, 3, 50, 150] do
    for m in [(2500000 : ℚ), 3300000, 4000000, 5000000, 9000000] do
      out := out ++ [(s!"nm/lambda={lam}/m={m}", { base with iso := isoD 1300 lam m })]
  return out

def oneOff : List (String × Inputs ℚ) :=
  [("base_rejected", base), ("base_pass", pass0),
   ("eps/wide", { base with eps := ⟨1/5, 1/10, 0, full1⟩ }),
   ("eps/sys_low", { base with eps := ⟨1/20, 1/50, 0, full1⟩ }),
   ("eps/omega_low", { base with eps := ⟨1/10, 1/100, 0, full1⟩ }),
   ("eps/omega_below_qspread", { base with eps := ⟨1/10, 1/50, 1/20, full1⟩ }),
   ("eps/omega_ge_1", { base with eps := ⟨1/10, 1, 0, full1⟩ }),
   ("val/imag_50", { base with val := { val0 with imagAbsMaxK := 50 } }),
   ("val/imag_signed_-500", { base with val := { val0 with imagAbsMaxK := -500 } }),
   ("val/HLiBe_no_check", { base with val := { val0 with containsHLiBe := true } }),
   ("val/HLiBe_check_shift_0.05", { base with
      val := { val0 with containsHLiBe := true, anharmonicCheck := true,
                         anharmShiftAbs := 1/20 },
      eps := ⟨1/10, 1/20, 0, full1⟩ }),
   ("val/shift_0.3_not_recomputed", { base with
      val := { val0 with containsHLiBe := true, anharmonicCheck := true,
                         anharmShiftAbs := 3/10 },
      eps := ⟨1/10, 3/10, 0, full1⟩ }),
   ("val/shift_0.3_recomputed", { base with
      val := { val0 with containsHLiBe := true, anharmonicCheck := true,
                         anharmShiftAbs := 3/10, alpha2FRecomputedAnharmonic := true },
      eps := ⟨1/10, 3/10, 0, full1⟩ }),
   ("val/shift_signed_-0.3", { base with
      val := { val0 with containsHLiBe := true, anharmonicCheck := true,
                         anharmShiftAbs := -3/10 } }),
   ("val/mu_negative", { base with val := { val0 with muStarNegativeSuspected := true } }),
   ("val/smearings_2", { base with val := { val0 with admittedSmearings := 2 } }),
   ("val/qmeshes_1", { base with val := { val0 with qMeshes := 1 } }),
   ("val/EF_missing", { base with val := { val0 with EF := none } }),
   ("val/EF_nonadiabatic", { base with val := { val0 with EF := some 5000 } }),
   ("val/EF_negative", { base with val := { val0 with EF := some (-5000) } }),
   ("val/nsheets_0", { base with val := { val0 with nSheets := 0 } }),
   ("val/nsheets_2_no_mb", { base with val := { val0 with nSheets := 2 } }),
   ("dom/A_zero", { base with iso := isoD 0 1 3000000 }),
   ("dom/A_negative", { base with iso := isoD (-1) 1 3000000 }),
   ("dom/lambda_negative", { base with iso := isoD 1250 (-1) 3000000 }),
   ("dom/dA_0.6", { base with iso := { (isoD 1250 1 3000000) with dA := 3/5 } }),
   ("dom/dA_negative", { base with iso := { (isoD 1250 1 3000000) with dA := -1/10 } }),
   ("dom/omega_max_negative", { base with val := { val0 with omegaMax := -1000 } }),
   ("dom/w_sum_0.9", xWsum),
   ("dom/w_lower_sum_gt_1", xWlow),
   ("dom/moments_inconsistent", { base with iso := isoD 1250 1 1000000 }),
   ("dom/w_n_mismatch", { base with dos := some (dos2 full1) }),
   ("dom/mat_B0", { base with mat := some ⟨0, fun i => i.elim0, fun i => i.elim0, fun i => i.elim0, 1/20, full1⟩ }),
   ("dom/mat_asymmetric", { base with mat := some ⟨2, ![![600, 100], ![300, 600]], ![1, 1], ![1, 1], 1/20, full1⟩ }),
   ("dom/mat_negative_entry", { base with mat := some ⟨2, ![![600, -10], ![-10, 600]], ![1, 1], ![1, 1], 1/20, full1⟩ }),
   ("dom/mat_v_zero", { base with mat := some ⟨2, ![![600, 100], ![100, 600]], ![1, 1], ![1, 0], 1/20, full1⟩ }),
   ("prov/eps_ML", { base with eps := ⟨1/10, 1/50, 0, pf .ml 1⟩ }),
   ("prov/eps_run2", { base with eps := ⟨1/10, 1/50, 0, pf .fullDFPT 2⟩ }),
   ("prov/lambda_coarse_m_external", xLamCoarse),
   ("prov/nmG_fires_full", { base with iso := isoD 1150 3 2000000 }),
   ("prov/A_external", { base with iso := isoD 500 1 1000000 .external }),
   ("tag/fileIncomplete", { base with tags := [.fileIncomplete] }),
   ("tag/fileCompleteUnverified", { base with tags := [.fileCompleteUnverified] }),
   ("tag/dropNonpos_0", { base with tags := [.dropNonpos 0] }),
   ("tag/dropNonpos_5", { base with tags := [.dropNonpos 5] }),
   ("tag/negClip_5", { base with tags := [.negClip 5] }),
   ("tag/sigmaExcluded_rej", { base with tags := [.sigmaExcluded true false] }),
   ("tag/sigmaExcluded_cert", { base with tags := [.sigmaExcluded false true] }),
   ("tag/informational", { base with tags := [.informational "M-ISO-SINGLEBAND"] }),
   ("tag/unknown", { base with tags := [.unknown "A-FILE-INCOMPLET"] }),
   ("tag/unknown_on_pass", { pass0 with tags := [.unknown "X"] }),
   ("kres/margin_pass_with_matrix", xKres),
   ("kres/margin_fail_no_matrix", { base with iso := isoD 1050 1 9000000 })]

def allCases : List (String × Inputs ℚ) := oneOff ++ lamCases ++ gridCases

def render : String :=
  let cases := allCases.map fun (id, x) => jcase cfg id x
  "{\n" ++
  "  \"generator\": \"lean/ExportCases.lean (TargetLogic.decideT at K = Q, exact rationals)\",\n" ++
  "  \"spec\": \"docs/theory/01_eliashberg_targets.md sec. 0.4 + 0.6 (post-RL1)\",\n" ++
  "  \"number_format\": \"exact rationals as strings p or p/q\",\n" ++
  "  \"config\": " ++ jcfg cfg ++ ",\n" ++
  "  \"n_cases\": " ++ toString cases.length ++ ",\n" ++
  "  \"cases\": [\n    " ++ ",\n    ".intercalate cases ++ "\n  ]\n}\n"

end ExportCases

#eval do
  IO.FS.createDirAll "export"
  IO.FS.writeFile "export/decide_cases.json" ExportCases.render
  IO.println s!"wrote {ExportCases.allCases.length} cases"
