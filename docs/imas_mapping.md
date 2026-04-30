# BPSD ↔ IMAS マッピングドラフト

BPSD の内部データ構造 (`bpsd_types.f90`) と IMAS Data Dictionary (IDS) の
対応関係を、移植・converter 設計のために洗い出すワーキング文書。
このドキュメントは「埋まる箇所」と「議論が必要な箇所」を区別することが目的で、
最終仕様ではない。

参考:

- IMAS-Python: <https://github.com/iterorganization/IMAS-Python>
- IMAS Data Dictionary (公開版): <https://imas-data-dictionary.readthedocs.io/>
- BPSD types: `bpsd_types.f90`

## 凡例

- ✅ 直接マッピング可（単位・座標も一致）
- ⚠️ 単位 / 座標 / 符号 / shape の調整が必要
- ❓ 開いている疑問（実データで突き合わせないと確定しない）
- ❌ IMAS 側に対応スロットが見当たらない / 概念が異なる

## 全体で要確認の事項

| 項目 | 内容 |
|---|---|
| `rho` の定義 | BPSD `%rho` は "normalized minor radius" だが rho_tor_norm (= sqrt(phi/phi_b)) か sqrt(psi_norm) か rho_pol か未確定。converter で IMAS の `grid.rho_tor_norm` / `grid.psi` のどれにマップするか決める必要あり |
| COCOS | IMAS は COCOS=11 推奨。BPSD は内部規約が明文化されていない（要 ats-fukuyama に確認）。Ip / Bt / psi 系の符号で問題化する |
| 時間 | BPSD: 各 `_type` に `%time` スカラー。IMAS: 各 IDS が `time(:)` 配列 + `time_slice(:)` のスライス。BPSD → IMAS は「単一スライスで書き出す」converter になる |
| species 識別 | BPSD: `bpsd_species_data` (pa, pz, npa) の配列で nsmax 種を区別。IMAS: `core_profiles.profiles_1d.electrons` と `ions(:)/element(:)/(a, z_n)` で電子と各イオンを別構造。converter は ns=1 (electron 規約?) を electrons に振り分ける必要あり |
| `idum` | 多くの type に `idum` がある（パディング用）。IMAS には対応概念なし、無視可 |

## bpsd_shot_type

`%shotID`, `%modelID`, `%deviceID` はメタデータで、IDS の中身ではなく
**IMAS Data Entry のキー** (pulse, run, machine) に対応する。

| BPSD | IMAS | メモ |
|---|---|---|
| `%shotID` | DataEntry の `pulse` | ✅ そのまま入る |
| `%modelID` | DataEntry の `run` | ⚠️ 概念が「モデル番号」と「ラン番号」でズレる場合あり |
| `%deviceID` | DataEntry の `machine` (または `dataset_description.ids_properties.provider`) | ✅ |

## bpsd_device_type

| BPSD | 単位 | IMAS path | IMAS 単位 | shape | 状態 |
|---|---|---|---|---|---|
| `%rr` (R0) | m | `equilibrium.vacuum_toroidal_field.r0` | m | scalar | ✅ |
| `%zz` | m | `equilibrium.time_slice/global_quantities/magnetic_axis.z` | m | per time | ⚠️ BPSD は「Geometrical vertical position」、IMAS は磁気軸 Z（同じか要確認） |
| `%ra` (a) | m | `equilibrium.time_slice/boundary/minor_radius` | m | per time | ✅ |
| `%rb` (wall) | m | `wall.description_2d(:)/limiter/unit(:)/outline.r/.z` | m | array | ❌ BPSD は典型値スカラー、IMAS は壁形状ポリゴン。「典型 wall radius」を IMAS に再現するスロットは無い |
| `%bb` (B0) | T | `equilibrium.vacuum_toroidal_field.b0(:)` | T | per time | ⚠️ 符号 (COCOS) |
| `%ip` | A | `equilibrium.time_slice/global_quantities/ip` | A | per time | ⚠️ 符号 (COCOS) |
| `%elip` | – | `equilibrium.time_slice/boundary/elongation` | – | per time | ✅ |
| `%trig` | – | `equilibrium.time_slice/boundary/triangularity_upper` & `triangularity_lower` | – | per time | ⚠️ BPSD は単一値、IMAS は upper/lower 別。対応は「平均値」か「両方に同じ値」か要決定 |

## bpsd_species_type / bpsd_species_data

| BPSD | IMAS | メモ |
|---|---|---|
| `%nsmax` | `core_profiles.profiles_1d/ions` の長さ + 1 (electrons 分) | ⚠️ ns=1 を電子と仮定するのか、ロジック不在 |
| `%data(ns)%pa` | `core_profiles.profiles_1d/ions(i)/element(:)/a` | ✅ 単位は同じ atomic mass number |
| `%data(ns)%pz` | `core_profiles.profiles_1d/ions(i)/z_ion` | ✅ |
| `%data(ns)%npa` | `core_profiles.profiles_1d/ions(i)/element(:)/z_n` | ✅ |

❓ 同位体 (D/T 混合) のように 1 ion species が複数 element からなるケースを
BPSD はサポートしていない（element 配列がない）。IMAS → BPSD では平均化、
BPSD → IMAS では element 配列長 1 で書き出す、になる。

## bpsd_equ1D_type

座標: BPSD `%rho`、IMAS `equilibrium.time_slice/profiles_1d/grid` 配下
(`rho_tor_norm`, `psi`, `phi`, `area`, `volume`)。

| BPSD | 単位 | IMAS path (`equilibrium.time_slice.profiles_1d/...`) | 状態 |
|---|---|---|---|
| `%data(:)%psit` | Wb | `phi` (toroidal flux) | ⚠️ BPSD 名 `psit` (toroidal) は IMAS の `phi` に対応。命名差注意。符号は COCOS |
| `%data(:)%psip` | Wb | `psi` (poloidal flux) | ⚠️ 同上、命名差 + COCOS |
| `%data(:)%ppp` | Pa | `pressure` | ✅ |
| `%data(:)%piq` (1/q) | – | `q` の逆数 | ⚠️ IMAS は `q` を持つ。converter で `q = 1.0/piq` (q→0 ハンドリング必要) |
| `%data(:)%pip` | A | `f` (= R·B_t) ? または積分系 | ❓ BPSD コメント "Poloidal current ~2πR·B/μ0" は実質 R·B_t。IMAS の `f` (= R·B_t、単位 T·m) と一致するか要確認、単位も違う |
| `%data(:)%pit` | A | `j_tor` の積分 / `equilibrium.time_slice/global_quantities/ip` | ❓ "Toroidal current ~2πr·Bp/μ0" がプロファイルとして何を意味するか要確認。j_tor profile を積分で再構築するのが自然 |

## bpsd_metric1D_type

IMAS `equilibrium.time_slice.profiles_1d` には flux-surface averaged metric
量が `gm1`～`gm9` として定義されている。BPSD の 20 個と機械的にマップ
できる訳ではない。

| BPSD | 単位 | IMAS path | 状態 |
|---|---|---|---|
| `%data(:)%pvol` | m^3 | `volume` | ✅ |
| `%data(:)%psur` | m^2 | `surface` | ✅ |
| `%data(:)%dvpsit` | m^3/Wb | `dvolume_dphi` | ✅ |
| `%data(:)%dvpsip` | m^3/Wb | `dvolume_dpsi` | ✅ |
| `%data(:)%aver2` (`<R^2>`) | m^2 | `gm5` または直接対応無し | ❓ IMAS の gm# 定義を実数で確認 |
| `%data(:)%aver2i` (`<1/R^2>`) | 1/m^2 | `gm1` (`<1/R^2>`) | ✅ 候補 |
| `%data(:)%aveb2` (`<B^2>`) | T^2 | `gm5` (`<B^2>`) | ✅ 候補 |
| `%data(:)%aveb2i` (`<1/B^2>`) | 1/T^2 | `gm4` (`<1/B^2>`) | ✅ 候補 |
| `%data(:)%avegv` (`<\|∇V\|>`) | m^2 | – | ❓ 直接スロット見当たらず |
| `%data(:)%avegv2` (`<\|∇V\|^2>`) | m^4 | `gm7` 系? | ❓ |
| `%data(:)%avegvr2` (`<\|∇V\|^2/R^2>`) | m^2 | `gm6`? | ❓ |
| `%data(:)%avegr` (`<\|∇ρ\|>`) | – | `gm7` | ✅ 候補 |
| `%data(:)%avegr2` (`<\|∇ρ\|^2>`) | – | `gm3` | ✅ 候補 |
| `%data(:)%avegrr2` (`<\|∇ρ\|^2/R^2>`) | 1/m^2 | `gm2` | ✅ 候補 |
| `%data(:)%avegpp2` (`<\|∇ψ_p\|^2>`) | Wb^2/m^2 | `gm8`? | ❓ |
| `%data(:)%rr` | m | `geometric_axis.r` | ✅ |
| `%data(:)%rs` | m | `r_outboard` - `r_inboard` の半分? | ⚠️ BPSD の "r" の正確な定義を要確認 |
| `%data(:)%elip` | – | `elongation` | ✅ |
| `%data(:)%trig` | – | `triangularity_upper/lower` | ⚠️ device と同じく upper/lower 分離 |
| `%data(:)%aveb` (`<B>`) | T | – | ❓ IMAS には `b_field_average` 系の直接スロット無し（gm# のいずれか） |

⚠️ gm1..gm9 の物理定義は IMAS DD 内の各 field の comment に記載がある。
converter を書く前に「BPSD の `<X>`」と「IMAS の `gmN`」の式が一致するか
**式単位で** 突合せる必要がある。

## bpsd_plasmaf_type

IMAS の対応先は `core_profiles.profiles_1d`。BPSD は `data(nr, ns)` の
2 次元配列で電子・各イオンを区別している。

ns 規約 (BPSD で ns=1 が電子か否か) は明文化されていない → ❓ 要確認。
以下では「ns=1=electron, ns≥2=ion(ns-1)」と仮定して書く。

| BPSD | 単位 | IMAS path (electron) | IMAS path (ion i) | 状態 |
|---|---|---|---|---|
| `%data(:,1)%density` | m^-3 | `electrons.density` | – | ✅ |
| `%data(:,i)%density` (i≥2) | m^-3 | – | `ions(i-1)/density` | ✅ |
| `%temperature` | eV | `electrons.temperature` | `ions(i-1)/temperature` | ⚠️ IMAS は通常 eV だが一部 IDS は J 系。要確認 |
| `%temperature_para` / `%temperature_perp` | eV | `electrons.temperature_fit?` 直接スロット限定的 | `ions(i-1)/.../temperature_para` 系 | ⚠️ IMAS の構造体の中で各方向が露出するパスを確認 |
| `%velocity_tor` | m/s | `electrons.velocity.toroidal` | `ions(i-1)/velocity.toroidal` | ✅ |
| `%velocity_pol` | m/s | `.../velocity.poloidal` | 同上 | ✅ |
| `%velocity_para` / `%velocity_perp` | m/s | `.../velocity.parallel` 系 | 同上 | ⚠️ |
| `%zave` | – | `core_profiles.profiles_1d.zeff`? | – | ❓ BPSD `zave` は「averaged charge」で species ごとに保持 (= z_ion?)。IMAS の zeff (= Σ n_i Z_i^2 / n_e) とは概念が違う。位置を要再検討 |
| `%z2ave` | – | – | – | ❌ IMAS に対応スロット見当たらず（converter 側で派生量として再計算） |
| `%density_fastion` | m^-3 | `core_profiles.profiles_1d.ions(i)/density_fast` | – | ✅ |
| `%energy_fastion` | eV | `core_profiles.profiles_1d.ions(i)/temperature_fast` 系 | – | ⚠️ 「fast ion energy」が average kinetic energy か total なのか要明確化 |
| `%qinv(nr)` | – | `equilibrium.time_slice.profiles_1d.q` の逆数 | – | ⚠️ IMAS では `q` は equilibrium 側にあって core_profiles 側ではない。BPSD では plasmaf に同梱 → 出力時にどちらに置くか方針決め |

## bpsd_trmatrix_type

IMAS 対応: `core_transport.model(:)/profiles_1d/...`

各 model は identifier (anomalous, neoclassical, ...) で分類される。
BPSD には model 区別が無く 1 セットのみ → converter では
「BPSD の trmatrix を anomalous + neoclassical のどちらに入れるか」あるいは
「unspecified model として 1 個書き出す」かの方針決定が必要 (❓)。

| BPSD | 単位 | IMAS path (`core_transport.model[X]/profiles_1d/electrons or ions(i)/...`) | 状態 |
|---|---|---|---|
| `%data(:,ns)%Dn` | m^2/s | `particles.d` | ✅ |
| `%data(:,ns)%un` | m/s | `particles.v` | ⚠️ IMAS では convection velocity の符号規約に注意 (内向きが正/負) |
| `%data(:,ns)%DT` | m^2/s | `energy.d` | ✅ |
| `%data(:,ns)%uT` | m/s | `energy.v` | ⚠️ 同上 |
| `%data(:,ns)%Dp` | m^2/s | `momentum_tor.d` ? | ❓ IMAS の momentum transport は方向別 (toroidal / parallel) で BPSD のスカラー `Dp` とそのまま噛み合わない |
| `%data(:,ns)%up` | m/s | `momentum_tor.v` ? | ❓ 同上 |

## bpsd_trsource_type

IMAS 対応: `core_sources.source(:)/profiles_1d/electrons or ions(i)/...`
各 source は `identifier` (NBI, ec_heating, ic_heating, lh_heating,
bremsstrahlung, line_radiation, cyclotron_radiation, ohmic, ionisation,
recombination, charge_exchange, ...) で分類される。

| BPSD | 単位 | IMAS source identifier | IMAS profiles path | 状態 |
|---|---|---|---|---|
| `%data(:,ns)%nip` | 1/(m^3 s) | `ionisation` | `electrons.particles` (e との対) / `ions(i)/particles` | ⚠️ ns との対応 (どの ion から?) を要決定 |
| `%data(:,ns)%nim` | 1/(m^3 s) | `recombination` | 同上 | ⚠️ 符号: BPSD は loss を正値で持つ可能性あり、IMAS は source 統一で sink は負 |
| `%data(:,ns)%ncx` | 1/(m^3 s) | `charge_exchange` | 同上 | ⚠️ 符号 |
| `%data(:,ns)%Pec` | W/m^3 | `ec` (electron cyclotron) | `electrons.energy` | ✅ |
| `%data(:,ns)%Plh` | W/m^3 | `lh` | `electrons.energy` | ✅ |
| `%data(:,ns)%Pic` | W/m^3 | `ic` | `ions(i)/energy` (主に) | ⚠️ ns で受け側を区別する規約を converter 側で明記 |
| `%data(:,ns)%Pbr` | W/m^3 | `bremsstrahlung` | `electrons.energy` | ⚠️ loss なので IMAS では負値で書き出すことになる (符号規約) |
| `%data(:,ns)%Pcy` | W/m^3 | `synchrotron_radiation` (旧 cyclotron radiation) | `electrons.energy` | ⚠️ 符号 |
| `%data(:,ns)%Plr` | W/m^3 | `line_radiation` | `electrons.energy` | ⚠️ 符号 |
| `%data(:,ns)%Poh` | W/m^3 | `ohmic` | `electrons.energy` | ✅ |

⚠️ 重要: BPSD の trsource は単位が `W1/m^3` と書かれている箇所があり
(古い trmatrix の流用コメント)、実際の単位を要確認。コード上は
`bpsd_setup_trsource_kdata` で `W/m^3` (W1 ではなく)。

## 未対応 / 後回し

| BPSD type | 状態 |
|---|---|
| `bpsd_equ2D_type` (psip 2D) | IMAS `equilibrium.time_slice.ggd` または `profiles_2d` に対応するが、現 BPSD では put/get 実装が薄い。後回し |
| `bpsd_equ3D_type` (psip 3D) | ステラレータ向け。IMAS では `equilibrium.time_slice.ggd` 系。後回し |
| `bpsd_dielectric_type` (3×3 complex tensor) | IMAS には `waves` や `core_profiles` の `velocity` などに分散して入る形で、直接スロットは無い。converter は対象外でよさそう |

## converter 設計に向けた次ステップ

1. **`%rho` の物理的定義の確定** (rho_tor_norm か sqrt(psi_norm) か)。
   ats-fukuyama 側で 1 件の実データを書き出して、IMAS の
   `equilibrium.time_slice.profiles_1d.rho_tor_norm` と数値比較するのが
   一番速い。
2. **BPSD の COCOS 規約の確認**。`%bb`, `%ip`, `%psip`, `%psit` の符号を
   1 ケースで実機データと突き合わせ、IMAS COCOS=11 との変換係数を確定。
3. **metric1D の gm# マッピングを式単位で検証**。IMAS DD の comment に
   ある積分式と BPSD の `aver2`, `avegv2`, ... の定義式を並べて確認する。
4. **trmatrix の momentum 取り扱い決定**。BPSD のスカラー Dp/up を IMAS
   の方向別 momentum transport にどう射影するか。
5. **prototype converter** を Python で実装 (BPSD save file `bpsd.data`
   または将来 NetCDF → IMAS-Python 経由で IDS に書き込む)。
   実装は最小から: `bpsd_device_type` → `equilibrium.vacuum_toroidal_field`
   と `bpsd_plasmaf_type` → `core_profiles.profiles_1d` の 2 経路だけで
   round-trip 1 件作る。
