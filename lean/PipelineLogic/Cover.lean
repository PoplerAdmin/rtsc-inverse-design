/-
P4: cover lemma for contiguous cells.
-/
import PipelineLogic.Common

namespace PipelineLogic

/-- `Tiles lo hi cells`: the list of cells `[(a₀,b₀), …, (a_K,b_K)]` is nonempty, `a₀ = lo`,
`b_i = a_{i+1}` (exact contiguity), and `b_K = hi`. (No ordering `a_i ≤ b_i` is needed.) -/
def Tiles : ℝ → ℝ → List (ℝ × ℝ) → Prop
  | _, _, [] => False
  | lo, hi, [c] => c.1 = lo ∧ c.2 = hi
  | lo, hi, c :: d :: cs => c.1 = lo ∧ Tiles c.2 hi (d :: cs)

/-- **P4 (cover lemma).** If `cells` tile `[lo, hi]` contiguously, every `T ∈ [lo, hi]` lies in
some cell `[a, b]`. -/
theorem cover : ∀ (cells : List (ℝ × ℝ)) (lo hi T : ℝ), Tiles lo hi cells →
    lo ≤ T → T ≤ hi → ∃ c ∈ cells, c.1 ≤ T ∧ T ≤ c.2
  | [], _, _, _, h, _, _ => h.elim
  | [c], _, _, _, ⟨h1, h2⟩, hl, hh => ⟨c, List.mem_singleton_self c, h1 ▸ hl, h2 ▸ hh⟩
  | c :: d :: cs, _, hi, T, ⟨h1, h2⟩, hl, hh => by
    by_cases hT : T ≤ c.2
    · exact ⟨c, List.mem_cons_self, h1 ▸ hl, hT⟩
    · obtain ⟨c', hc', h3, h4⟩ := cover (d :: cs) c.2 hi T h2 (le_of_lt (not_le.mp hT)) hh
      exact ⟨c', List.mem_cons_of_mem _ hc', h3, h4⟩

/-- **P4 (point-list form).** Breakpoints `t 0 = lo, …, t K = hi`: every `T ∈ [lo, hi]` lies
in some `[t i, t (i+1)]` with `i < K`, provided `K ≥ 1`. -/
theorem cover_points (t : ℕ → ℝ) : ∀ (K : ℕ), 1 ≤ K → ∀ T, t 0 ≤ T → T ≤ t K →
    ∃ i < K, t i ≤ T ∧ T ≤ t (i + 1)
  | 0, h, _, _, _ => absurd h (by norm_num)
  | K + 1, _, T, h0, hK => by
    by_cases hT : t K ≤ T
    · exact ⟨K, Nat.lt_succ_self K, hT, hK⟩
    · have hK1 : 1 ≤ K := by
        rcases Nat.eq_zero_or_pos K with h | h
        · subst h; exact absurd h0 hT
        · exact h
      obtain ⟨i, hi, h1, h2⟩ := cover_points t K hK1 T h0 (le_of_lt (not_le.mp hT))
      exact ⟨i, Nat.lt_succ_of_lt hi, h1, h2⟩

end PipelineLogic
