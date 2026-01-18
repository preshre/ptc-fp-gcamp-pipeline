"""
PTC Fiber Photometry Analysis - Main Script

Signal: 465nm (GCaMP) - Calcium-dependent fluorescence
Control: 405nm (Isobestic) - Motion artifact correction

Outputs: dF/F, PETH, AUC, Ymax
"""

import os
import argparse
import logging
from datetime import datetime

from pl_fp_mouse_465 import Mouse
from utils.utils import (
    resolve_long_path, 
    parse_text_file_auc, 
    convert_to_excel_with_formatting_auc,
    parse_text_file_ymax, 
    convert_to_excel_with_formatting_ymax
)


def get_timestamp():
    """Get current timestamp for results folder naming."""
    return datetime.now().strftime("%Y%m%d_%H%M")


def process_single_mouse(block_path, signal, control, result_dir, is_train, 
                         pre_time, post_time, auc_log_path, ymax_log_path):
    """
    Process a single mouse's fiber photometry data.
    
    Parameters:
    -----------
    block_path : str
        Path to TDT data block
    signal : str
        Signal channel name (465nm GCaMP)
    control : str
        Control channel name (405nm isobestic)
    result_dir : str
        Output directory for results
    is_train : bool
        True for training, False for LTM
    pre_time : int
        Pre-stimulus window (seconds)
    post_time : int
        Post-stimulus window (seconds)
    auc_log_path : str
        Path for AUC log file
    ymax_log_path : str
        Path for Ymax log file
    """
    mouse_name = os.path.basename(block_path)
    print(f'\n{"="*50}')
    print(f'Processing: {mouse_name}')
    print(f'Signal: {signal} (465nm GCaMP)')
    print(f'Control: {control} (405nm Isobestic)')
    print(f'{"="*50}')
    
    mouse = Mouse(
        file_path=block_path,
        is_train=is_train,
        pre_time=pre_time,
        post_time=post_time,
        result_dir=result_dir,
        signal=signal,
        control=control,
        auc_log_path=auc_log_path,
        ymax_log_path=ymax_log_path
    )
    
    mouse.run_analysis(mouse_name)


def process_directory(dir_path, signal="_465A", control="_405A"):
    """
    Process all TDT data blocks in a directory.
    
    Parameters:
    -----------
    dir_path : str
        Directory containing TDT data
    signal : str
        Signal channel (default: "_465A" for 465nm GCaMP)
    control : str
        Control channel (default: "_405A" for 405nm isobestic)
    """
    # Create results directory
    timestamp = get_timestamp()
    results_dir = os.path.join("Results", f"PTC_Results_{timestamp}")
    os.makedirs(results_dir, exist_ok=True)
    
    # Set up log file paths
    auc_log_path = resolve_long_path(os.path.join(results_dir, 'auc.txt'))
    auc_xlsx_path = resolve_long_path(os.path.join(results_dir, 'auc.xlsx'))
    ymax_log_path = resolve_long_path(os.path.join(results_dir, 'yMax.txt'))
    ymax_xlsx_path = resolve_long_path(os.path.join(results_dir, 'yMax.xlsx'))
    
    # Remove existing log files
    for path in [auc_log_path, ymax_log_path]:
        if os.path.exists(path):
            os.remove(path)
    
    print(f'\n{"#"*60}')
    print(f'PTC Fiber Photometry Analysis')
    print(f'Signal: 465nm (GCaMP) | Control: 405nm (Isobestic)')
    print(f'Outputs: dF/F, PETH, AUC, Ymax')
    print(f'Results directory: {results_dir}')
    print(f'{"#"*60}\n')
    
    # Process all TDT blocks in directory
    processed_count = 0
    for root, _, files in os.walk(dir_path):
        root_lower = root.lower()
        
        # Skip archive and habituation folders
        if any(keyword in root_lower for keyword in ('archive', 'habituation')):
            logging.info(f"Skipping: {root}")
            continue
        
        # Check for TDT files (.tev)
        if any(file.endswith('.tev') for file in files):
            is_train = 'train' in root_lower
            
            try:
                process_single_mouse(
                    block_path=root,
                    signal=signal,
                    control=control,
                    result_dir=results_dir,
                    is_train=is_train,
                    pre_time=30,
                    post_time=60.1,
                    auc_log_path=str(auc_log_path),
                    ymax_log_path=str(ymax_log_path)
                )
                processed_count += 1
            except Exception as e:
                logging.error(f"Error processing {root}: {e}")
                print(f"Error processing {root}: {e}")
    
    # Generate Excel output files
    print(f'\n{"="*50}')
    print('Generating Excel output files...')
    
    try:
        if os.path.exists(auc_log_path):
            data = parse_text_file_auc(str(auc_log_path))
            convert_to_excel_with_formatting_auc(data, str(auc_xlsx_path))
            print(f'AUC data saved to: {auc_xlsx_path}')
    except Exception as e:
        print(f"Error generating AUC Excel: {e}")
    
    try:
        if os.path.exists(ymax_log_path):
            data = parse_text_file_ymax(str(ymax_log_path))
            convert_to_excel_with_formatting_ymax(data, str(ymax_xlsx_path))
            print(f'Ymax data saved to: {ymax_xlsx_path}')
    except Exception as e:
        print(f"Error generating Ymax Excel: {e}")
    
    print(f'\n{"#"*60}')
    print(f'Analysis complete! Processed {processed_count} recordings.')
    print(f'Results saved to: {results_dir}')
    print(f'{"#"*60}\n')


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="PTC Fiber Photometry Analysis\n"
                    "Signal: 465nm (GCaMP) | Control: 405nm (Isobestic)\n"
                    "Outputs: dF/F, PETH, AUC, Ymax"
    )
    parser.add_argument(
        "--signal", 
        type=str, 
        default="_465A", 
        help="Signal channel name (default: '_465A' for 465nm GCaMP)"
    )
    parser.add_argument(
        "--control", 
        type=str, 
        default="_405A", 
        help="Control channel name (default: '_405A' for 405nm isobestic)"
    )
    parser.add_argument(
        "--dir",
        type=str,
        default=None,
        help="Directory path containing TDT data (optional, will prompt if not provided)"
    )
    
    args = parser.parse_args()
    
    # Get directory path
    if args.dir:
        dir_path = args.dir
    else:
        dir_path = input("Enter the directory path containing TDT data: ").strip()
    
    dir_path = resolve_long_path(dir_path)
    
    if os.path.isdir(dir_path):
        process_directory(str(dir_path), signal=args.signal, control=args.control)
    else:
        print("Error: The provided path is not a valid directory.")


if __name__ == "__main__":
    main()
