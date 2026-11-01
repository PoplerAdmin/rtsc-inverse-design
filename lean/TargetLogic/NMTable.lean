/-
TargetLogic.NMTable: the certified NM-G table (doc 01 §0.3, [THM-C]).
GENERATED from docs/theory/certificates/nm_*.json (exact decimal `m0_over_T2`; do not edit).
Entry `(λ_k, m0_k)`: certificate Theorem NM(λ_k, m0_k) (CERT-NM); g_lo(λ_k) = 90000·m0_k at 300 K.
Source files and sha256:
  nm_0.8.json  ac983e23bac6c5d9d5e215235550d98035c3e39a33c90f2faae67206782aa456
  nm_0.9.json  3d06eb5592aa0acd98022c904f6397e5a72eff317e89e0723ec4b98cae44161d
  nm_1.0.json  ea03029f17074ebf1d300ad5cdf8cc6982a110e17f154669e8802968803a58d8
  nm_1.1.json  07699a0ea7d9f864dd6071803f9929f631b7c22437151e4faa0d1dcf39b17ad3
  nm_1.2.json  b6f24824c811ade1eaf4cce009cde752ab1c2722eae5effc338ce6733120a36b
  nm_1.3.json  393226620c45d4df3614fd2c0b38783ec336b9f71cc5bbde5f122d134c6b1412
  nm_1.4.json  ed424913e4befb343019ce8e41a8b22c7387a5cc54b13659f6f19dcf9f32c3a2
  nm_1.5.json  917cb26c91d67f9ec18991ae2188852953fbc9a7a0824a11f28d01d061c7cc42
  nm_1.6.json  71a28ac649c1e96438f10b9ec3c6ead0eaa6c6866d2cb9ad84f64be5ede9724a
  nm_1.7.json  b3466e2fe95cbcaf50ca390bbd296b0a7f50728786691958b90b95a17d6fe64d
  nm_1.8.json  38617208d6b549d45310b32673f61364e19a5efe61c3e3d311605878223d9f2d
  nm_1.9.json  085ee30213d51ae2a69e0f169e53d455c3967bc26ffeca293d1bc665c671f7d6
  nm_2.0.json  22eb38c5f567169d61c799adce76261ed150f498a2ff7f9e4cc9c4b39232be80
  nm_2.25.json  6bdf794e6baac32946e5979e3b233e64225dea6de62410b29732573cecb4ab7c
  nm_2.5.json  a92b5f055acdb330a47b249c5a762bd2f594744517ef5b3dd01092629a8da7bc
  nm_2.75.json  d9cb24b6f29de5f2b758391a1547e8d073744b5d2e661014cba16b41a81eefcb
  nm_3.0.json  2010bd339a6e93df8b14836d8816cc17bd8399c69b6d92912d1ef93ca6505e35
  nm_3.25.json  4b4fbf9b4ce28f0c60c04f3f6da53c7c1462927004ef56506f8726344a41f1f8
  nm_3.5.json  931cbb98729f9642e333cf8b14a6558d6c392dad0dccdeff2b25a4256ee2b470
  nm_3.75.json  f7263df341d5b3a59e5f12c63ba6f66aa64cbbc518627326f5dcbc29e2fbe039
  nm_4.0.json  4639a65ab1dcaf1bc93ba2cb26ec19a2cf5bb5f22909fb012cfd54aa0bbcda92
  nm_4.5.json  35375ceba53d294b3756016bf9bf8dacf20042c7163dc4418e855bdb9e2c21da
  nm_5.0.json  380eaee471b864a9ef15d2608a188febb793e158bc80e073e76a6af182376247
  nm_6.0.json  3ffa3fc837fcd037fe394d24a17584947ade4f19628451c440868760490b0ee9
  nm_7.0.json  7401494eb37ddc9b8b08c3714ebea9de624b49ae0bf98a39bf9525defcd0b684
  nm_8.0.json  86d0bc5f0602eed5a32c4446982b6747b564d81fff37844e53bafd52bffb2340
  nm_10.0.json  972ff4c58a9b466b696e5ab79c4240f430e53d7363d1e45f1595e115de91ec38
  nm_15.0.json  481629c51e8e01727db22f8668d6c88a9c2e0f723d1f4a30770816b4316eeb3c
  nm_20.0.json  eee6e3705af0728224ae7b39693b3addebe79a8616a709783db3dcb76f4f76e9
  nm_50.0.json  b75e94a4668e2150fae71e8c8ec1aaf5302e51f88ab0b18c181001398c728f4d
  nm_100.0.json  997dae775bacc9b004cb2fdd20417cc88bbb10c677a046c1ff74c45bfeb2c9d6
  nm_inf.json  7e92cae144a715b5aafde1a87f78e2991135a023f914b16f1219c6cfa7fe85d7
-/
import Mathlib.Algebra.Order.Field.Basic

namespace TargetLogic

section
variable {K : Type*} [Field K]

/-- `(λ_k, m0_k)` from the certificates, sorted by λ_k. -/
def nmM0Table : List (K × K) :=
  [((4 / 5 : K), (1050297903 / 10000000 : K)),
   ((9 / 10 : K), (8746099489 / 100000000 : K)),
   ((1 / 1 : K), (1522329739 / 20000000 : K)),
   ((11 / 10 : K), (3415051477 / 50000000 : K)),
   ((6 / 5 : K), (391511279 / 6250000 : K)),
   ((13 / 10 : K), (46704013 / 800000 : K)),
   ((7 / 5 : K), (5506799827 / 100000000 : K)),
   ((3 / 2 : K), (327669819 / 6250000 : K)),
   ((8 / 5 : K), (5027630707 / 100000000 : K)),
   ((17 / 10 : K), (2424650727 / 50000000 : K)),
   ((9 / 5 : K), (58739947 / 1250000 : K)),
   ((19 / 10 : K), (4571198223 / 100000000 : K)),
   ((2 / 1 : K), (4460820829 / 100000000 : K)),
   ((9 / 4 : K), (2120900177 / 50000000 : K)),
   ((5 / 2 : K), (4079233901 / 100000000 : K)),
   ((11 / 4 : K), (1976955351 / 50000000 : K)),
   ((3 / 1 : K), (1927204547 / 50000000 : K)),
   ((13 / 4 : K), (235845411 / 6250000 : K)),
   ((7 / 2 : K), (3706502513 / 100000000 : K)),
   ((15 / 4 : K), (1825033521 / 50000000 : K)),
   ((4 / 1 : K), (720380269 / 20000000 : K)),
   ((9 / 2 : K), (1762025873 / 50000000 : K)),
   ((5 / 1 : K), (3463870771 / 100000000 : K)),
   ((6 / 1 : K), (337692837 / 10000000 : K)),
   ((7 / 1 : K), (165858031 / 5000000 : K)),
   ((8 / 1 : K), (3273556149 / 100000000 : K)),
   ((10 / 1 : K), (642840327 / 20000000 : K)),
   ((15 / 1 : K), (784500091 / 25000000 : K)),
   ((20 / 1 : K), (1550552853 / 50000000 : K)),
   ((50 / 1 : K), (1518315461 / 50000000 : K)),
   ((100 / 1 : K), (3015674861 / 100000000 : K))]

/-- `m0` of the λ = ∞ certificate. -/
def nmM0Inf : K := (2994979021 / 100000000 : K)

end

end TargetLogic
