from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Fragments
import math

def count_hydroxylic(mol):

    return Fragments.fr_Al_OH(mol)

def count_carboxylic(mol):

    return Fragments.fr_Al_COO(mol)


def count_primary_secondary_amines(mol):
    """第一級および第二級アミンの数を数える"""

    return Fragments.fr_NH2(mol) + Fragments.fr_NH1(mol)

def count_amide(mol):
    """アミド結合 (-CONH-) の数を数える"""
    return Fragments.fr_amide(mol)

def count_sulfoxide_bond(mol):
    """S=O 結合 (スルホキシドなど) の数を数える"""
    # 元のコードは S=O をカウントしていたため、関数名を修正
    count = 0
    pattern = Chem.MolFromSmarts('[#16]=[#8]') # S=O の SMARTS パターン
    if pattern:
        matches = mol.GetSubstructMatches(pattern)
        count = len(matches)
    return count


def count_ring(mol):
    """縮合していない環と縮合環の数を数える"""
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings() # 環を構成する原子のインデックスのリストのリスト
    num_rings = len(atom_rings)

    # RDKit の GetAtomRings は最小の環のセット (SSSR) とは限らないため、
    # 縮合の判定が複雑になることがある。
    # ここでは、より簡単なアプローチとして、
    # 環を共有する原子を持つ環を縮合環とみなす。
    # ただし、スピロ結合などは単純には扱えない。

    # より堅牢な方法：環を構成する結合を調べる
    bond_rings = ring_info.BondRings()
    fused_rings_count = 0
    all_ring_atoms = set()
    for ring in atom_rings:
        all_ring_atoms.update(ring)

    # 縮合に関与する原子（複数の環に属する原子）を数える
    fused_atoms = set()
    if num_rings > 1:
         for i in range(num_rings):
             for j in range(i + 1, num_rings):
                 shared_atoms = set(atom_rings[i]).intersection(set(atom_rings[j]))
                 if len(shared_atoms) > 0: # 原子を共有していれば縮合の可能性
                    fused_atoms.update(shared_atoms)
                    # 厳密には結合を共有しているかで判断すべき
                    # ここでは簡略化のため原子共有で判断

    # RDKitの GetRingInfo は SSSR ではないため、単純な数え方では
    # 複雑な縮合系で問題が出る可能性がある。
    # 例えばナフタレンは2つの環として数えられるが、
    # より複雑な分子では意図しない数え方になることも。

    # 元のコードのロジックを踏襲しつつ、原子ベースで判定
    # (ただし、元のコードは BondRings を使っており、より正確だった可能性がある)
    # 元のコードのロジック（修正版）：
    fused_rings_indices = set()
    for i, ring_bonds in enumerate(bond_rings):
        is_fused = False
        for bond_idx in ring_bonds:
            bond = mol.GetBondWithIdx(bond_idx)
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            #結合の両端の原子が、それぞれ複数の環に属しているかチェック
            if ring_info.NumAtomRings(begin_atom.GetIdx()) > 1 and \
               ring_info.NumAtomRings(end_atom.GetIdx()) > 1:
                 # 共有結合が複数の環に属しているかで判定する方がより正確
                 if ring_info.NumBondRings(bond_idx) > 1:
                      is_fused = True
                      break
        if is_fused:
             fused_rings_indices.add(i)

    num_fused_rings = len(fused_rings_indices)
    num_unfused_rings = num_rings - num_fused_rings

    # 安全のため、マイナスにならないようにする
    num_unfused_rings = max(0, num_unfused_rings)

    return num_unfused_rings, num_fused_rings


def count_frag(mol):
    """分子内の特徴的なフラグメントを数える"""
    num_hydroxylic = count_hydroxylic(mol)
    num_carboxylic = count_carboxylic(mol)
    num_primary_secondary_amines = count_primary_secondary_amines(mol) # 関数名を修正
    num_amide = count_amide(mol)
    num_sulfoxide_bond = count_sulfoxide_bond(mol) # 関数名を修正
    num_unfused_rings, num_fused_rings = count_ring(mol)
    return {
        'num_hydroxylic': num_hydroxylic,
        'num_carboxylic': num_carboxylic,
        'num_primary_secondary_amines': num_primary_secondary_amines,
        'num_amide': num_amide,
        'num_sulfoxide_bond': num_sulfoxide_bond, # キー名を修正
        'num_unfused_rings': num_unfused_rings,
        'num_fused_rings': num_fused_rings
    }

def estimate_density(mol):
    """
    分子オブジェクトから密度 (g/cm^3) を推定する。
    元のコードのロジックを基にする。
    """
    if mol is None:
        return None # 計算不可

    frag = count_frag(mol)
    # 元のコードの補正項。係数は経験的なものと思われる。
    correction = (frag['num_hydroxylic'] * 0.1 +
                  frag['num_carboxylic'] * 0.1 +
                  frag['num_primary_secondary_amines'] * 0.1 +
                  frag['num_amide'] * 0.1 +
                  frag['num_sulfoxide_bond'] * 0.1 + # キー名を修正
                  frag['num_unfused_rings'] * 0.1 +
                  frag['num_fused_rings'] * 0.075)

    # 補正項の上限
    correction = min(correction, 0.3)

    # 水素原子を付加して分子量と Vs を計算
    mol_with_hs = Chem.AddHs(mol)
    if mol_with_hs is None:
        # 水素付加に失敗した場合 (サニタイズ問題など)
        return None

    # 分子量 M (g/mol)
    M = sum([a.GetMass() for a in mol_with_hs.GetAtoms()])

    # Vs (原子価殻電子に基づく経験的な体積項？) の計算
    Vs = 0
    atomic_numbers = [a.GetAtomicNum() for a in mol_with_hs.GetAtoms()]
    for n in atomic_numbers:
        if n <= 2:        # Period 1 (H, He)
            Vs += 1
        elif n >= 3 and n <= 10:  # Period 2 (Li to Ne)
            Vs += 2
        elif n >= 11 and n <= 18: # Period 3 (Na to Ar)
            Vs += 4
        elif n >= 19 and n <= 36: # Period 4 (K to Kr)
            Vs += 5
        elif n >= 37 and n <= 54: # Period 5 (Rb to Xe)
            Vs += 7.5
        elif n >= 55 and n <= 86: # Period 6 (Cs to Rn)
            Vs += 9
        else:             # Period 7 onwards
            Vs += 9 # 元コードに合わせてPeriod 6と同じ値を使用

    # 元の密度計算式。係数 5.0 も経験的な値と思われる。
    # M/Vs はおおよそ (質量/体積のような項) に対応
    if Vs == 0: return None # ゼロ除算回避
    density_g_cm3 = (1 + correction) * M / Vs / 5.0

    return density_g_cm3

# --- 新しい関数 ---
def calculate_cell_length(smiles: str, num_molecules: int):
    """
    SMILES文字列とセル内の分子数を指定して、
    推定密度に基づき立方体セルの辺長 (nm) を計算する。

    Args:
        smiles (str): 分子のSMILES文字列。
        num_molecules (int): シミュレーションセルに入れる分子の数。

    Returns:
        float: 推定されたセルの一辺の長さ (nm)。
               計算できなかった場合は None を返す。
    Raises:
        ValueError: 無効なSMILES文字列または分子数が0以下の場合。
    """
    if num_molecules <= 0:
        raise ValueError("分子数は正の整数である必要があります。")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"無効なSMILES文字列です: {smiles}")

    # 分子の密度を推定 (g/cm^3)
    density_g_cm3 = estimate_density(mol)
    if density_g_cm3 is None or density_g_cm3 <= 0:
         # 密度が計算できないか、非物理的な値の場合はエラーまたはNoneを返す
         print(f"警告: {smiles} の密度を推定できませんでした。")
         return None # またはエラーを発生させる

    # 分子量 (g/mol) - 再度計算が必要 (estimate_density内で計算しているが返り値ではないため)
    mol_with_hs = Chem.AddHs(mol)
    if mol_with_hs is None:
        print(f"警告: {smiles} に水素を付加できませんでした。")
        return None
    M = sum([a.GetMass() for a in mol_with_hs.GetAtoms()])

    # アボガドロ定数 (mol^-1)
    N_A = 6.02214076e23

    density_g_cm3 = density_g_cm3 - 0.1

    # セルの体積 V (cm^3) を計算
    # V = (N * M) / (rho * N_A)
    volume_cm3 = (num_molecules * M) / (density_g_cm3 * N_A)

    # セルの一辺の長さ L (cm) を計算
    length_cm = volume_cm3**(1/3)

    # セル長を nm に変換 (1 cm = 1e7 nm)
    length_nm = length_cm * 1e7

    return length_nm, density_g_cm3


if __name__ == '__main__':
    import sys
    import argparse # 引数処理をより分かりやすくするために argparse を使用

    parser = argparse.ArgumentParser(description='Estimate cubic cell length (nm) for a given SMILES and number of molecules.')
    parser.add_argument('smiles', type=str, help='SMILES string of the molecule.')
    parser.add_argument('num_molecules', type=int, help='Number of molecules in the cell.')

    if len(sys.argv) == 1:
        # 引数がない場合はヘルプを表示（あるいはエラーでも良い）
        parser.print_help(sys.stderr)
        sys.exit(1)
        
    args = parser.parse_args()

    try:
        # calculate_cell_length 関数を呼び出し
        cell_length_nm, density_g_cm3  = calculate_cell_length(args.smiles, args.num_molecules)
        #print(f"Estimate_density(mol): {density_g_cm3 + 0.1 :.6f} g/cm^3", file=sys.stderr)

        if cell_length_nm is not None:
            # 成功したら結果を標準出力に表示
            print(f"{cell_length_nm:.6f}") # 十分な精度で出力
            sys.exit(0) # 正常終了
        else:
            # 関数が None を返した場合 (内部でエラー処理済みのはず)
            print(f"Error: Could not calculate cell length for SMILES '{args.smiles}'. Density estimation might have failed.", file=sys.stderr)
            sys.exit(1) # 異常終了

    except ValueError as e:
        # calculate_cell_length が ValueError を投げた場合 (無効な SMILES など)
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1) # 異常終了
    except Exception as e:
        # その他の予期せぬエラー
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1) # 異常終了