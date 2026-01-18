"""
PTC Fiber Photometry Analysis - Mouse Class

Signal: 465nm (GCaMP) - Calcium-dependent fluorescence
Control: 405nm (Isobestic) - Motion artifact correction

Outputs: dF/F, PETH, AUC, Ymax
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt, savgol_filter
from scipy.optimize import curve_fit
from scipy.integrate import trapezoid
from scipy.stats import sem
import tdt
import os
import math
import warnings
import seaborn as sns

# Suppress overflow warnings during curve fitting (expected behavior with large time values)
warnings.filterwarnings("ignore")

from utils.extract_timestamps import extract_timestamps_from_excel


class Mouse:
    """
    Class for fiber photometry data analysis of individual mice.
    
    Performs signal processing and generates dF/F, PETH, AUC, and Ymax outputs.
    """
    
    def __init__(self, 
                 file_path: str, 
                 is_train: bool, 
                 pre_time: int, 
                 post_time: int, 
                 result_dir: str,
                 signal: str = "_465A", 
                 control: str = "_405A", 
                 auc_log_path: str = 'auc.txt',
                 ymax_log_path: str = 'yMax.txt'):
        """
        Initialize Mouse object and process fiber photometry data.
        
        Parameters:
        -----------
        file_path : str
            Path to TDT data block
        is_train : bool
            True for training session, False for LTM retrieval
        pre_time : int
            Pre-stimulus window duration in seconds
        post_time : int
            Post-stimulus window duration in seconds
        result_dir : str
            Directory for saving output files
        signal : str
            Signal channel name (default: "_465A" for 465nm GCaMP)
        control : str
            Control channel name (default: "_405A" for 405nm isobestic)
        auc_log_path : str
            Path for AUC log file
        ymax_log_path : str
            Path for Ymax log file
        """
        self.BLOCK_PATH = file_path
        self.PRE_TIME = pre_time
        self.POST_TIME = post_time
        self.is_train = is_train
        self.auc_log_path = auc_log_path
        self.ymax_log_path = ymax_log_path
        self.result_dir = result_dir
        self.is_err = False
        
        # Extract mouse name from path
        self.mouse_name = self.BLOCK_PATH.split('\\')[-1].split('_')[0].split('-')[0]
        if '/' in self.BLOCK_PATH:
            self.mouse_name = self.BLOCK_PATH.split('/')[-1].split('_')[0].split('-')[0]
        
        self.experiment_name = self._determine_experiment_name(self.BLOCK_PATH)
        
        # Initialize data attributes
        self.t1 = None
        self.data = None
        self.time = None
        self.fs = None
        self.CSon = None
        self.CSoff = None
        self.USon = None
        self.USoff = None
        self.signal_lowpass = None
        self.control_lowpass = None
        self.dFF = None
        self.zScore = None
        self.savgol_zscore = None
        self.dFF_snips = None
        self.peri_time = None
        
        # Timestamp storage
        self.training_timestamps = None
        self.ltm1d_timestamps = None
        self.ltm14d_timestamps = None
        self.ltm28d_timestamps = None
        
        # Load and process data
        print(f'Loading data for mouse: {self.mouse_name}')
        self._load_data()
        self._calculate_properties(signal, control)
    
    def _load_data(self):
        """Load TDT fiber photometry data."""
        data = tdt.read_block(self.BLOCK_PATH)
        
        def load_timestamps():
            """Load timestamps from Excel files."""
            script_dir = os.path.dirname(os.path.abspath(__file__))
            default_file = os.path.join(script_dir, "PTC Timestamps_CS_US.xlsx")
            
            timestamps = extract_timestamps_from_excel(default_file)
            (self.training_timestamps, self.ltm1d_timestamps, 
             self.ltm14d_timestamps, self.ltm28d_timestamps, _, _) = timestamps
            
            # Check for mouse-specific timestamp file
            mouse_file_upper = os.path.join(script_dir, f"PTC Timestamps_CS_US {self.mouse_name.upper()}.xlsx")
            mouse_file_lower = os.path.join(script_dir, f"PTC Timestamps_CS_US {self.mouse_name.lower()}.xlsx")
            
            if os.path.isfile(mouse_file_upper):
                ts = extract_timestamps_from_excel(mouse_file_upper)
                self._update_timestamps(ts)
            elif os.path.isfile(mouse_file_lower):
                ts = extract_timestamps_from_excel(mouse_file_lower)
                self._update_timestamps(ts)
        
        def get_experiment_start_time():
            """Get start time based on experiment type."""
            if self.experiment_name == "LTM28d":
                return self.ltm28d_timestamps["Start"][0]
            elif self.experiment_name == "LTM14d":
                return self.ltm14d_timestamps["Start"][0]
            return self.ltm1d_timestamps["Start"][0]
        
        # Check for valid epocs data
        if hasattr(data.epocs, 'PrtB') and hasattr(data.epocs, 'CS__'):
            self.t1 = data.epocs.PrtB.onset[0]
            t2 = data.epocs.CS__.offset[1] + 1800 if self.is_train else 0
            self.data = tdt.read_block(self.BLOCK_PATH, t1=self.t1, t2=t2)
            return
        
        # Handle error case - use timestamps from Excel
        self.is_err = True
        load_timestamps()
        
        if self.is_train:
            self.t1 = self.training_timestamps["Start"][0]
            t2 = self.training_timestamps["CSoff"][-1] + 1800
        else:
            self.t1 = get_experiment_start_time()
            t2 = 0
        
        self.data = tdt.read_block(self.BLOCK_PATH, t1=self.t1, t2=t2)
    
    def _update_timestamps(self, timestamps):
        """Update timestamps if valid."""
        training, ltm1d, ltm14d, ltm28d, _, _ = timestamps
        
        def validate(new_ts, default_ts):
            if new_ts is not None and not math.isnan(new_ts['CSon'][0]):
                return new_ts
            return default_ts
        
        self.training_timestamps = validate(training, self.training_timestamps)
        self.ltm1d_timestamps = validate(ltm1d, self.ltm1d_timestamps)
        self.ltm14d_timestamps = validate(ltm14d, self.ltm14d_timestamps)
        self.ltm28d_timestamps = validate(ltm28d, self.ltm28d_timestamps)
    
    def _calculate_properties(self, signal, control):
        """
        Process fiber photometry signals and calculate derived metrics.
        
        Signal Processing Pipeline:
        1. Low-pass filtering (6Hz Butterworth)
        2. Photobleaching correction (double exponential fit)
        3. Motion correction (linear fit using control channel)
        4. dF/F calculation
        5. Z-score calculation with Savitzky-Golay smoothing
        """
        # Get signal and control data
        self.signal = self.data.streams[signal].data
        self.control = self.data.streams[control].data
        
        # Ensure equal lengths
        max_len = min(len(self.signal), len(self.control))
        self.signal = self.signal[:max_len]
        self.control = self.control[:max_len]
        
        # Calculate sampling rate and time vector
        self.fs = self.data.streams[signal].fs
        self.time = np.linspace(1, len(self.signal), len(self.signal)) / self.fs
        
        # Extract CS/US timestamps
        self._extract_timestamps()
        
        # Step 1: Low-pass filtering
        self.signal_lowpass, self.control_lowpass = self._lowpass_filter(
            self.signal, self.control, self.fs
        )
        
        # Step 2: Photobleaching correction
        signal_corrected, control_corrected, signal_fitted, _ = self._double_fit_subtraction(
            self.signal_lowpass, self.control_lowpass, self.time
        )
        
        # Step 3: Motion correction and dF/F calculation
        self.dFF, Y_dF_all, _ = self._linear_fit_subtraction(
            signal_corrected, control_corrected, signal_fitted
        )
        
        # Step 4: Z-score calculation
        self.zScore = (Y_dF_all - np.mean(Y_dF_all)) / np.std(Y_dF_all)
        self.savgol_zscore = self._get_savgol_zscore(5000)
        
        # Step 5: Calculate peri-event snippets
        self.dFF_snips, self.peri_time = self._calculate_dFF_snips(
            self.savgol_zscore, self.fs, self.CSon, self.time, 
            self.PRE_TIME, self.POST_TIME
        )
        
        # Log Ymax values
        ymax_values = [round(np.max(snippet), 2) for snippet in self.dFF_snips]
        if self.is_train:
            US_snips, _ = self._calculate_dFF_snips(
                self.savgol_zscore, self.fs, self.CSon, self.time, -28, 30.1
            )
            ymax_values.extend([round(np.max(snippet), 2) for snippet in US_snips])
        self._log_ymax(ymax_values)
    
    def _extract_timestamps(self):
        """Extract CS and US onset/offset times."""
        if self.is_err:
            if self.is_train:
                self.USon = np.array(self.training_timestamps['USon']) - self.t1
                self.USoff = np.array(self.training_timestamps['USoff']) - self.t1
                self.CSon = np.array(self.training_timestamps['CSon']) - self.t1
                self.CSoff = np.array(self.training_timestamps['CSoff']) - self.t1
            else:
                timestamps = {
                    "LTM28d": self.ltm28d_timestamps,
                    "LTM14d": self.ltm14d_timestamps,
                }.get(self.experiment_name, self.ltm1d_timestamps)
                self.CSon = np.array(timestamps["CSon"]) - self.t1
                self.CSoff = np.array(timestamps["CSoff"]) - self.t1
        else:
            if self.is_train:
                self.USon = self.data['epocs']['Shck']['onset'] - self.t1
                self.USoff = self.data['epocs']['Shck']['offset'] - self.t1
            self.CSon = self.data['epocs']['CS__']['onset'] - self.t1
            self.CSoff = self.data['epocs']['CS__']['offset'] - self.t1
    
    def _lowpass_filter(self, signal, control, fs):
        """
        Apply 6Hz low-pass Butterworth filter.
        
        Parameters:
        -----------
        signal : array
            465nm GCaMP signal
        control : array
            405nm isobestic control
        fs : float
            Sampling frequency
            
        Returns:
        --------
        tuple : Filtered signal and control arrays
        """
        order = 3
        cutoff = 6
        nyquist = fs / 2
        normalized_cutoff = cutoff / nyquist
        
        b, a = butter(order, normalized_cutoff, btype='low')
        filtered_signal = filtfilt(b, a, signal)
        filtered_control = filtfilt(b, a, control)
        
        return filtered_signal, filtered_control
    
    def _double_fit_subtraction(self, signal, control, time):
        """
        Perform photobleaching correction using double exponential fit.
        
        Fits: y = a * exp(b*x) + c * exp(d*x)
        
        Returns:
        --------
        tuple : Corrected signal, corrected control, fitted signal curve, fitted control curve
        """
        def double_exponential(x, a, b, c, d):
            return a * np.exp(b * x) + c * np.exp(d * x)
        
        signal = np.array(signal, dtype=np.float32)
        control = np.array(control, dtype=np.float32)
        time = np.array(time, dtype=np.float32)
        
        popt_signal, _ = curve_fit(double_exponential, time, signal, 
                                   p0=[0.1, -0.1, 0.1, -0.1], maxfev=10000)
        popt_control, _ = curve_fit(double_exponential, time, control, 
                                    p0=[0.1, -0.1, 0.1, -0.1], maxfev=10000)
        
        fitted_signal = double_exponential(time, *popt_signal)
        fitted_control = double_exponential(time, *popt_control)
        
        signal_corrected = signal - fitted_signal
        control_corrected = control - fitted_control
        
        return signal_corrected, control_corrected, fitted_signal, fitted_control
    
    def _linear_fit_subtraction(self, signal, control, fitted_signal):
        """
        Perform motion correction using linear fit.
        
        Uses 405nm isobestic control to remove motion artifacts from 465nm signal.
        
        Returns:
        --------
        tuple : dF/F values, delta F values, fitted values
        """
        bls = np.polyfit(control, signal, 1)
        fitted_values = bls[0] * control + bls[1]
        delta_F = signal - fitted_values
        
        epsilon = 1e-9
        dFF = 100 * delta_F / (fitted_signal + epsilon)
        
        return dFF, delta_F, fitted_values
    
    def _get_savgol_zscore(self, window_length=5000):
        """Apply Savitzky-Golay smoothing to z-score."""
        return savgol_filter(self.zScore, window_length, polyorder=3)
    
    def _calculate_dFF_snips(self, zscore, fs, CSon, time, pre_time, post_time):
        """
        Extract peri-event time histogram snippets.
        
        Parameters:
        -----------
        zscore : array
            Z-scored signal
        fs : float
            Sampling frequency
        CSon : array
            CS onset times
        time : array
            Time vector
        pre_time : int
            Pre-stimulus window
        post_time : int
            Post-stimulus window
            
        Returns:
        --------
        tuple : Array of snippets, peri-event time vector
        """
        trange = np.array([-pre_time, post_time]) * fs
        trange = trange.astype(int)
        snippet_length = trange[1] - trange[0]
        
        dFF_snips = []
        for on in CSon:
            if on >= pre_time:
                onset_idx = np.searchsorted(time, on)
                pre_idx = onset_idx + trange[0]
                post_idx = onset_idx + trange[1]
                
                if pre_idx >= 0 and post_idx < len(zscore):
                    dFF_snips.append(zscore[pre_idx:post_idx])
                else:
                    dFF_snips.append(np.zeros(snippet_length))
        
        peri_time = np.linspace(-pre_time, post_time, snippet_length)
        return np.array(dFF_snips), peri_time
    
    def _determine_experiment_name(self, filename):
        """Determine experiment type from filename."""
        filename_lower = filename.lower()
        if 'ltm28d' in filename_lower or 'ltm_28d' in filename_lower:
            return 'LTM28d'
        elif 'ltm14d' in filename_lower or 'ltm_14d' in filename_lower:
            return 'LTM14d'
        elif 'ltm1' in filename_lower or 'ltm_1d' in filename_lower:
            return 'LTM1d'
        return 'Training'
    
    def _log_auc(self, auc_values):
        """Log AUC values to file."""
        with open(self.auc_log_path, 'a') as f:
            f.write(f'{self.mouse_name}-{self.experiment_name}: {auc_values}\n')
    
    def _log_ymax(self, ymax_values):
        """Log Ymax values to file."""
        with open(self.ymax_log_path, 'a') as f:
            f.write(f'{self.mouse_name}-{self.experiment_name}: {ymax_values}\n')
    
    def _save_plot(self, fig, output_dir, filename):
        """Save plot as PNG and PDF."""
        animal_id = filename.split('_')[0].split('-')[0]
        experiment_name = self._determine_experiment_name(filename)
        
        output_dir_lower = output_dir.lower()
        if 'peth' in output_dir_lower:
            file_suffix = 'PETH'
        elif 'auc' in output_dir_lower:
            file_suffix = 'AUC'
        elif 'dff' in output_dir_lower:
            file_suffix = 'dFF'
        else:
            file_suffix = 'plot'
        
        if self.is_train:
            base_filename = f'{animal_id}_train_{file_suffix}'
        else:
            base_filename = f'{animal_id}_{experiment_name}_{file_suffix}'
        
        output_path = os.path.join(self.result_dir, output_dir, animal_id)
        os.makedirs(output_path, exist_ok=True)
        
        fig.savefig(os.path.join(output_path, f'{base_filename}.png'), dpi=300, bbox_inches='tight')
        fig.savefig(os.path.join(output_path, f'{base_filename}.pdf'), dpi=300, bbox_inches='tight')
        plt.close(fig)
    
    def plot_dFF(self, filename):
        """
        Plot dF/F trace.
        
        Output: dF/F trace showing normalized fluorescence change over time.
        """
        # Determine max time for plotting
        max_time = 570 if self.is_train else 720
        if not self.is_err:
            try:
                end_time = self.data.epocs.End_.onset[-1] - self.t1
                if end_time > 570:
                    max_time = float(end_time)
            except:
                pass
        
        indices = np.where(self.time > max_time)[0]
        last_index = indices[0] if indices.size > 0 else len(self.time)
        
        fig, ax = plt.subplots(figsize=(12, 4))
        ax.plot(self.time[:last_index], self.dFF[:last_index], 
                linewidth=1, color='black', alpha=0.7)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
        ax.tick_params(axis='both', which='major', labelsize=8)
        plt.rc('font', family='Arial')
        ax.set_ylabel(r"$\Delta$F/F")
        ax.set_xlabel('Time (s)')
        plt.margins(x=0, y=0)
        plt.tight_layout()
        
        self._save_plot(fig, 'dFF', filename)
        print('Saved dF/F plot')
    
    def plot_PETH(self, filename):
        """
        Plot Peri-Event Time Histogram.
        
        Output: Mean response with SEM aligned to CS onset.
        """
        mean_response = np.mean(self.dFF_snips, axis=0)
        sem_response = sem(self.dFF_snips, axis=0)
        
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # Plot mean with SEM
        ax.plot(self.peri_time, mean_response, linewidth=1, color='black', 
                label='Mean Response', rasterized=True)
        ax.fill_between(self.peri_time, mean_response + sem_response, 
                        mean_response - sem_response,
                        facecolor='black', alpha=0.3, label='SEM', rasterized=True)
        
        # Plot CS marker
        ymin = np.min(self.dFF_snips) - 0.5
        ymax = np.max(self.dFF_snips) + 0.5
        diff = ymax - ymin
        ax.set_ylim(ymin, ymax)
        ax.plot([0, 30], [ymax, ymax], color='#BB1F1F', linewidth=2, label='CS', rasterized=True)
        
        if self.is_train:
            ax.plot([28, 30], [ymax + (0.02 * diff), ymax + (0.02 * diff)], color='#040404', 
                    linewidth=2, label='US', rasterized=True)
        
        ax.set_xlabel('Time (s)', fontname='Arial', fontsize=10)
        ax.set_ylabel(r'$\Delta$F/F', fontname='Arial', fontsize=10)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_xticks(np.arange(min(self.peri_time), max(self.peri_time) + 1, 10))
        ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.1f'))
        ax.tick_params(axis='both', which='major', labelsize=8)
        plt.rc('font', family='Arial')
        plt.margins(x=0, y=0)
        plt.tight_layout()
        
        self._save_plot(fig, 'PETH Histograms', filename)
        print('Saved PETH plot')
    
    def calculate_and_plot_AUC(self, filename):
        """
        Calculate and plot Area Under Curve.
        
        Time intervals:
        - Pre-CS: -30 to 0 seconds
        - CS: 0 to 30 seconds  
        - Post-CS: 30 to 60 seconds
        
        Output: Bar plot with AUC values and SEM error bars.
        """
        time_intervals = [(-30, 0), (0, 30), (30, 60)]
        num_events = len(self.dFF_snips)
        
        auc_pre = np.zeros(num_events)
        auc_cs = np.zeros(num_events)
        auc_post = np.zeros(num_events)
        
        for i, dFF in enumerate(self.dFF_snips):
            for start, end in time_intervals:
                start_idx = np.where(self.peri_time >= start)[0][0]
                end_idx = np.where(self.peri_time >= end)[0][0]
                time_subset = self.peri_time[start_idx:end_idx + 1]
                auc = trapezoid(dFF[start_idx:end_idx + 1], time_subset)
                
                if start == -30:
                    auc_pre[i] = auc
                elif start == 0:
                    auc_cs[i] = auc
                else:
                    auc_post[i] = auc
        
        # Log AUC values
        logging_array = [
            np.round(auc_pre, 2).tolist(),
            np.round(auc_cs, 2).tolist(),
            np.round(auc_post, 2).tolist()
        ]
        self._log_auc(logging_array)
        
        # Create bar plot
        fig, ax = plt.subplots()
        x_labels = ['Pre-CS', 'CS', 'Post-CS']
        avg_values = [np.mean(auc_pre), np.mean(auc_cs), np.mean(auc_post)]
        sem_values = [sem(auc_pre), sem(auc_cs), sem(auc_post)]
        
        ax.bar(x_labels, avg_values, yerr=sem_values, capsize=5, 
               color="#BB1F1F", width=0.8)
        ax.set_xlabel("CS Event")
        ax.set_ylabel("AUC Values")
        ax.set_title("AUC for CS Events" if self.is_train else "AUC - LTM")
        plt.tight_layout()
        
        self._save_plot(fig, 'AUC Plots', filename)
        print('Saved AUC plot')
        
        return auc_pre, auc_cs, auc_post
    
    def plot_heat_map(self, size_x, size_y, filename, title="", ylabel=""):
        """
        Plot a heatmap of the delta F/F snippets using Seaborn.

        Plots a heatmap of the delta F/F snippets calculated for each CS event,
        with time on the x-axis and CS event on the y-axis.

        Parameters:
        - size_x (int): Width of the heatmap plot.
        - size_y (int): Height of the heatmap plot.
        - title (str): Title of the plot.
        - ylabel (str): Label for the y-axis (default is an empty string).
        """
        # Calculate actual time values for x-axis
        time_range = np.linspace(-self.PRE_TIME, self.POST_TIME, len(self.dFF_snips[0]))

        # Define x-tick labels at 10-second intervals
        xtick_labels = np.arange(-self.PRE_TIME, self.POST_TIME + 1, 10)

        # Create the figure and axis for the heatmap
        fig, ax = plt.subplots(figsize=(size_x, size_y))

        # Reverse the order of dFF_snips for display purposes
        dFF_snips_reversed = np.flip(self.dFF_snips, axis=0)

        # Generate the heatmap with Seaborn
        sns.heatmap(dFF_snips_reversed, cmap='jet', ax=ax, cbar=True)

        # Set labels and title
        ax.set_ylabel(ylabel)
        ax.set_title(title)

        # Map the x-ticks to the appropriate time points
        xtick_positions = np.interp(xtick_labels, time_range, np.arange(len(self.dFF_snips[0])))
        ax.set_xticks(xtick_positions)
        ax.set_xticklabels([f"{int(label)}" for label in xtick_labels])

        # Set the y-ticks and reverse the labels for correct CS event display
        ax.set_yticks(np.arange(len(self.dFF_snips)) + 0.5)
        ax.set_yticklabels(np.arange(1, len(self.dFF_snips) + 1)[::-1])

        # Set the axis limits
        ax.set_xlim(0, len(self.dFF_snips[0]) - 1)
        ax.set_ylim(0, len(self.dFF_snips))

        self._save_plot(fig, 'Heatmaps', filename)
        print('########## Saved Heatmap ##########')
    
    def run_analysis(self, filename):
        """
        Run complete analysis pipeline.
        
        Outputs:
        - dF/F trace plot
        - PETH plot with mean ± SEM
        - AUC bar plot
        - Logged Ymax values
        """
        print(f'\n*** Starting analysis for {self.mouse_name} ***')
        
        self.plot_dFF(filename)
        self.plot_PETH(filename)
        self.calculate_and_plot_AUC(filename)
        self.plot_heat_map(8, 4, filename, title="", ylabel="")
        
        print(f'*** Analysis complete for {self.mouse_name} ***\n')
