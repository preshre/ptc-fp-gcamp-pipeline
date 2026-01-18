"""
Extract CS/US timestamps from Excel files for PTC experiments.
"""

import pandas as pd


def extract_timestamps_from_excel(file_path):
    """
    Extract timestamps for Training, LTM1d, LTM14d, and LTM28d from Excel file.
    
    Parameters:
    -----------
    file_path : str
        Path to timestamps Excel file
        
    Returns:
    --------
    tuple : (training_timestamps, ltm1d_timestamps, ltm14d_timestamps, 
             ltm28d_timestamps, None, None)
    """
    def extract_timestamps(table, is_training=False):
        """Extract CS onset/offset times from a table."""
        timestamps = {"CSon": [], "CSoff": []}
        
        for _, row in table.iterrows():
            epoch = str(row[1])
            if epoch == 'Start':
                timestamps['Start'] = (row[2], row[3])
            elif epoch.startswith('CS+'):
                timestamps["CSon"].append(row[2])
                timestamps["CSoff"].append(row[3])
            elif epoch == 'End':
                timestamps['End'] = (row[2], row[3])
        
        # Calculate US times for training (US occurs 2s before CS offset)
        if is_training:
            timestamps["USon"] = [cs_off - 2 for cs_off in timestamps["CSoff"]]
            timestamps["USoff"] = timestamps["CSoff"]
        
        return timestamps
    
    # Load Excel file
    sheet_data = pd.read_excel(file_path, sheet_name=0, header=None)
    
    # Find table start indices
    training_start = sheet_data[sheet_data[0] == 'Training'].index[0]
    ltm1d_start = sheet_data[sheet_data[0] == 'LTM1d'].index[0]
    ltm14d_start = sheet_data[sheet_data[0] == 'LTM14d'].index[0]
    ltm28d_start = sheet_data[sheet_data[0] == 'LTM28d'].index[0]
    
    # Extract Training timestamps
    training_table = sheet_data.iloc[training_start + 2:ltm1d_start].dropna(how='all')
    training_table.columns = sheet_data.iloc[training_start + 1]
    training_table = training_table.reset_index(drop=True)
    training_timestamps = extract_timestamps(training_table, is_training=True)
    
    # Extract LTM1d timestamps
    ltm1d_table = sheet_data.iloc[ltm1d_start + 2:ltm14d_start].dropna(how='all')
    ltm1d_table.columns = sheet_data.iloc[ltm1d_start + 1]
    ltm1d_table = ltm1d_table.reset_index(drop=True)
    ltm1d_timestamps = extract_timestamps(ltm1d_table, is_training=False)
    
    # Extract LTM14d timestamps
    ltm14d_table = sheet_data.iloc[ltm14d_start + 2:ltm28d_start].dropna(how='all')
    ltm14d_table.columns = sheet_data.iloc[ltm14d_start + 1]
    ltm14d_table = ltm14d_table.reset_index(drop=True)
    ltm14d_timestamps = extract_timestamps(ltm14d_table, is_training=False)
    
    # Extract LTM28d timestamps
    ltm28d_table = sheet_data.iloc[ltm28d_start + 2:].dropna(how='all')
    ltm28d_table.columns = sheet_data.iloc[ltm28d_start + 1]
    ltm28d_table = ltm28d_table.reset_index(drop=True)
    ltm28d_timestamps = extract_timestamps(ltm28d_table, is_training=False)
    
    return (training_timestamps, ltm1d_timestamps, ltm14d_timestamps, 
            ltm28d_timestamps, None, None)
