"""
Memory Tracking and Logging Module for SpecMap

This module provides utilities to track and log memory usage during HSI data processing.
It helps identify memory bottlenecks and prevents out-of-memory crashes with large datasets.

Key features:
- Track peak memory usage per operation
- Log memory usage with GB-level granularity
- Provide warnings when memory thresholds are exceeded
- Always flush logs immediately to disk for crash diagnostics

Usage:
    from memory_tracker import MemoryTracker
    
    tracker = MemoryTracker()
    
    # Track an operation
    with tracker.track("Loading spectra"):
        # ... load data ...
        pass
    
    # Or manually
    tracker.log_memory("Before loading", context="Data loading")
    # ... do work ...
    tracker.log_memory("After loading", context="Data loading")
"""

import os
import gc
import psutil
import logging
from datetime import datetime
from contextlib import contextmanager


class MemoryTracker:
    """
    Track and log memory usage for SpecMap operations.
    
    This class provides utilities to monitor memory consumption during data processing
    and log it to a file for diagnostic purposes. Logs are flushed immediately after
    each write to ensure they are available even if the application crashes.
    
    Attributes:
        log_file (str): Path to the memory log file
        logger (logging.Logger): Logger instance for memory tracking
        process (psutil.Process): Process object for memory monitoring
        warning_threshold_gb (float): Memory threshold in GB for warnings
    """
    
    def __init__(self, log_file='logs/memory_usage.log', warning_threshold_gb=8.0):
        """
        Initialize the memory tracker.
        
        Args:
            log_file (str): Path to the memory log file. Defaults to 'logs/memory_usage.log'
            warning_threshold_gb (float): Threshold in GB to trigger warnings (default: 8.0)
        """
        self.log_file = log_file
        self.warning_threshold_gb = warning_threshold_gb
        self.process = psutil.Process()
        
        # Create logs directory if it doesn't exist
        log_dir = os.path.dirname(log_file)
        if log_dir and not os.path.exists(log_dir):
            try:
                os.makedirs(log_dir, exist_ok=True)
            except Exception as e:
                print(f"Warning: Could not create logs directory '{log_dir}': {e}")
                self.log_file = 'memory_usage.log'
        
        # Configure logger
        self.logger = logging.getLogger('MemoryTracker')
        self.logger.setLevel(logging.DEBUG)
        
        # Prevent duplicate handlers if re-initialized
        if not self.logger.handlers:
            # Create file handler with immediate flush
            try:
                file_handler = logging.FileHandler(
                    self.log_file,
                    mode='a',
                    encoding='utf-8'
                )
                file_handler.setLevel(logging.DEBUG)
                
                # Create formatter with timestamp and memory info
                formatter = logging.Formatter(
                    '[%(asctime)s] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S'
                )
                file_handler.setFormatter(formatter)
                
                # Add handler to logger
                self.logger.addHandler(file_handler)
                
            except Exception as e:
                print(f"Error setting up memory logging: {e}")
        
        # Log initialization
        self._log_header()
    
    def _log_header(self):
        """Log a header with system information."""
        total_ram = psutil.virtual_memory().total / (1024**3)
        self.logger.info("=" * 80)
        self.logger.info(f"Memory Tracking Session Started")
        self.logger.info(f"System Total RAM: {total_ram:.2f} GB")
        self.logger.info(f"Warning Threshold: {self.warning_threshold_gb:.2f} GB")
        self.logger.info("=" * 80)
        self._flush()
    
    def _flush(self):
        """Flush all handlers to ensure logs are written immediately."""
        for handler in self.logger.handlers:
            handler.flush()
    
    def get_memory_usage(self):
        """
        Get current memory usage statistics.
        
        Returns:
            dict: Dictionary with memory usage in bytes and GB
                - 'rss_bytes': Resident Set Size in bytes
                - 'rss_gb': Resident Set Size in GB
                - 'available_gb': Available system memory in GB
                - 'percent': Percentage of system memory used
        """
        # Force garbage collection to get accurate reading
        gc.collect()
        
        # Get process memory info
        mem_info = self.process.memory_info()
        
        # Get system memory info
        sys_mem = psutil.virtual_memory()
        
        return {
            'rss_bytes': mem_info.rss,
            'rss_gb': mem_info.rss / (1024**3),
            'available_gb': sys_mem.available / (1024**3),
            'percent': sys_mem.percent
        }
    
    def log_memory(self, operation, context="", data_info=None):
        """
        Log current memory usage for a specific operation.
        
        Args:
            operation (str): Name/description of the operation
            context (str): Additional context information
            data_info (dict): Optional data about what's being processed
                             (e.g., {'spectra_count': 1000, 'points_per_spectrum': 1024})
        
        Returns:
            dict: Current memory usage statistics
        """
        mem = self.get_memory_usage()
        
        # Build log message
        msg_parts = [
            f"[{operation}]",
            f"Memory: {mem['rss_gb']:.3f} GB",
            f"Available: {mem['available_gb']:.2f} GB",
            f"({mem['percent']:.1f}% system)"
        ]
        
        if context:
            msg_parts.insert(1, f"Context: {context}")
        
        if data_info:
            info_str = ", ".join([f"{k}={v}" for k, v in data_info.items()])
            msg_parts.append(f"| Data: {info_str}")
        
        log_message = " | ".join(msg_parts)
        
        # Check if we should warn
        if mem['rss_gb'] >= self.warning_threshold_gb:
            self.logger.warning(f"  HIGH MEMORY: {log_message}")
        else:
            self.logger.info(log_message)
        
        # Always flush immediately
        self._flush()
        
        return mem
    
    def log_memory_delta(self, operation, mem_before, context=""):
        """
        Log memory change since a previous measurement.
        
        Args:
            operation (str): Name of the operation
            mem_before (dict): Previous memory stats from get_memory_usage()
            context (str): Additional context
        
        Returns:
            dict: Current memory usage statistics
        """
        mem_after = self.get_memory_usage()
        delta_gb = mem_after['rss_gb'] - mem_before['rss_gb']
        
        # Build log message
        sign = "+" if delta_gb >= 0 else ""
        msg = (
            f"[{operation}] Memory: {mem_after['rss_gb']:.3f} GB "
            f"({sign}{delta_gb:.3f} GB) | Available: {mem_after['available_gb']:.2f} GB"
        )
        
        if context:
            msg += f" | Context: {context}"
        
        if delta_gb >= 1.0:
            self.logger.warning(f"  LARGE INCREASE: {msg}")
        elif mem_after['rss_gb'] >= self.warning_threshold_gb:
            self.logger.warning(f"  HIGH MEMORY: {msg}")
        else:
            self.logger.info(msg)
        
        # Always flush immediately
        self._flush()
        
        return mem_after
    
    @contextmanager
    def track(self, operation, context="", data_info=None):
        """
        Context manager to track memory usage for an operation.
        
        Args:
            operation (str): Name of the operation
            context (str): Additional context
            data_info (dict): Optional data information
        
        Example:
            with tracker.track("Loading spectra", data_info={'count': 1000}):
                # ... load data ...
                pass
        """
        # Log before
        mem_before = self.log_memory(f"{operation} - START", context=context, data_info=data_info)
        
        try:
            yield
        finally:
            # Log after
            self.log_memory_delta(f"{operation} - END", mem_before, context=context)
    
    def log_separator(self, title=""):
        """Log a separator line for readability."""
        if title:
            self.logger.info(f"{'='*40} {title} {'='*40}")
        else:
            self.logger.info("="*80)
        self._flush()


# Module-level convenience function
_default_memory_tracker = None

def get_default_memory_tracker():
    """
    Get or create the default MemoryTracker instance.
    
    This provides a singleton tracker that can be safely accessed
    from any thread or context in the application.
    
    Returns:
        MemoryTracker: The default memory tracker instance
    """
    global _default_memory_tracker
    if _default_memory_tracker is None:
        _default_memory_tracker = MemoryTracker()
    return _default_memory_tracker


# Example usage
if __name__ == "__main__":
    import numpy as np
    
    print("Testing MemoryTracker...")
    
    tracker = MemoryTracker(log_file='test_memory.log')
    
    # Test 1: Basic logging
    tracker.log_separator("Test 1: Basic Logging")
    tracker.log_memory("Initial state")
    
    # Test 2: Track with context manager
    tracker.log_separator("Test 2: Context Manager")
    with tracker.track("Creating large array", data_info={'size': '100MB'}):
        large_array = np.random.randn(10000, 1024)
    
    # Test 3: Manual delta tracking
    tracker.log_separator("Test 3: Manual Delta")
    mem_before = tracker.get_memory_usage()
    another_array = np.random.randn(5000, 1024)
    tracker.log_memory_delta("Created another array", mem_before)
    
    # Cleanup
    del large_array, another_array
    gc.collect()
    tracker.log_memory("After cleanup")
    
    print(f"\nMemory log written to: {tracker.log_file}")
