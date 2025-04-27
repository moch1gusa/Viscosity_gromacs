# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# --- 設定項目 ---
DEFAULT_INPUT_FILE = "viscosity_einstein.xvg"
DEFAULT_OUTPUT_FILE = "viscosity_avg_plot.png" # 出力ファイル名変更

# --- .xvg ファイル内の列インデックス ---
# head コマンドの結果、5列データと判明。最後の列が平均粘度と仮定。
TIME_COL_IDX = 0     # 時間 t の列インデックス
AVG_VISC_COL_IDX = -1 # 平均粘度 η_avg(t) の列インデックス (最後の列)
# ----------------------------------------

def plot_average_viscosity(input_filename, output_filename):
    """
    GROMACSの viscosity_einstein.xvg (5列データと仮定) を読み込み、
    最後の列 (平均粘度と仮定) vs 時間 のグラフをPNGで保存する関数
    Args:
        input_filename (str): 入力となる .xvg ファイルのパス
        output_filename (str): グラフの出力ファイルパス
    """
    print(f"入力ファイル: {input_filename}")
    print(f"出力ファイル: {output_filename}")

    # --- ファイル読み込み ---
    try:
        data = np.loadtxt(input_filename, comments=['@', '#'])
        print(f"データ読み込み完了。データ形状: {data.shape}")
    except FileNotFoundError:
        print(f"エラー: 入力ファイル '{input_filename}' が見つかりません。")
        sys.exit(1)
    except ValueError as e:
        print(f"エラー: ファイル '{input_filename}' のデータ形式が不正です。{e}")
        sys.exit(1)
    except Exception as e:
        print(f"エラー: ファイル '{input_filename}' の読み込み中に予期せぬエラー: {e}")
        sys.exit(1)

    # --- データ検証 ---
    actual_avg_visc_col_idx = AVG_VISC_COL_IDX if AVG_VISC_COL_IDX >= 0 else data.shape[1] + AVG_VISC_COL_IDX
    if data.ndim != 2 or data.shape[1] <= max(TIME_COL_IDX, actual_avg_visc_col_idx) or min(TIME_COL_IDX, actual_avg_visc_col_idx)<0:
        print(f"エラー: データの次元または列数が不足しています。形状: {data.shape}")
        print(f"       スクリプトは時間(列{TIME_COL_IDX})と平均粘度(列{actual_avg_visc_col_idx})を想定しています。")
        sys.exit(1)
    if data.shape[0] < 2:
        print("エラー: データ点が少なすぎます（2点未満）。")
        sys.exit(1)

    # --- データ抽出 ---
    try:
        time = data[:, TIME_COL_IDX]
        viscosity_avg = data[:, actual_avg_visc_col_idx] # 最後の列

        valid_indices = time >= 0
        if not np.all(valid_indices):
             print("警告: 時間が負のデータ点が見つかりました。これらは除外されます。")
             time = time[valid_indices]
             viscosity_avg = viscosity_avg[valid_indices]
             if time.size < 2:
                  print("エラー: 有効なデータ点が少なすぎます（2点未満）。")
                  sys.exit(1)
    except IndexError:
         print(f"エラー: 指定された列インデックス ({TIME_COL_IDX}, {actual_avg_visc_col_idx}) がデータ範囲外です。形状: {data.shape}")
         sys.exit(1)
    except Exception as e:
         print(f"エラー: データ抽出中に予期せぬエラーが発生しました: {e}")
         sys.exit(1)

    # --- グラフ描画 ---
    plt.figure(figsize=(10, 6)) # サイズ調整 (元の12,7でもOK)
    plt.plot(time, viscosity_avg, label=f'Average Viscosity (Column {actual_avg_visc_col_idx})')

    # --- グラフの装飾 ---
    plt.xlabel("Time (ps)")
    plt.ylabel("Average Viscosity ($Pa \cdot s$)")
    plt.title(f"Average Viscosity vs Time ({os.path.basename(input_filename)})")
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.6)

    current_ylim = plt.ylim()
    if current_ylim[0] < 0:
       plt.ylim(bottom=0)

    # --- グラフ保存 ---
    try:
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        print(f"グラフを '{output_filename}' として正常に保存しました。")
        plt.close()
    except Exception as e:
        print(f"エラー: グラフの保存中にエラーが発生しました: {e}")
        sys.exit(1)

# --- メイン処理 ---
if __name__ == "__main__":
    input_file = DEFAULT_INPUT_FILE
    output_file = "" # Initialize output_file

    if len(sys.argv) > 1:
        input_file = sys.argv[1]

    base, ext = os.path.splitext(input_file)
    # Use the new default suffix
    output_file = base + DEFAULT_OUTPUT_FILE[-len("_visc_plot.png"):] if DEFAULT_OUTPUT_FILE.endswith("_visc_plot.png") else base + "_avg_visc_plot.png"

    if len(sys.argv) > 2:
        output_file = sys.argv[2]

    plot_average_viscosity(input_file, output_file)
