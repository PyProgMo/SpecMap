"""
Centralized Error Handling and Logging Module for SpecMap Application

This module provides centralized error handling and logging functionality for the SpecMap application.
It includes a rotating file handler for persistent logging and user-friendly error messaging via tkinter.

Example Usage:
    # Initialize the error handler at application startup
    from error_handler import ErrorEngine
    
    error_engine = ErrorEngine()
    
    # Use in try/except blocks
    try:
        # some operation
        pass
    except Exception as e:
        error_engine.error(
            exception=e,
            context="Reading spectrum file",
            filename="example.txt"
        )
    
    # Get logger for direct logging
    logger = error_engine.get_logger()
    logger.info("Application started successfully")
"""

import logging
from logging.handlers import RotatingFileHandler
import os
import traceback
from datetime import datetime

# Try to import tkinter, but make it optional for testing
try:
    import tkinter as tk
    HAS_TKINTER = True
except ImportError:
    HAS_TKINTER = False
    # Create dummy messagebox for environments without tkinter
    class messagebox:
        @staticmethod
        def showerror(title, message):
            print(f"ERROR [{title}]: {message}")
        
        @staticmethod
        def showwarning(title, message):
            print(f"WARNING [{title}]: {message}")


class ContextFilter(logging.Filter):
    """Filter to ensure context attribute exists in all log records."""
    def filter(self, record):
        if not hasattr(record, 'context'):
            record.context = 'Unknown'
        return True


class ErrorEngine:
    """
    Centralized error handling engine for SpecMap application.
    
    Handles exceptions by:
    1. Logging detailed information to a rotating log file
    2. Displaying user-friendly error messages via tkinter messagebox
    3. Providing thread-safe logging capabilities
    
    Attributes:
        logger (logging.Logger): Configured logger instance
        log_file (str): Path to the log file
    """
    
    def __init__(self, log_file='logs/specmap.log', max_bytes=10*1024*1024, backup_count=5):
        """
        Initialize the ErrorEngine with rotating file handler.
        
        Args:
            log_file (str): Path to the log file. Defaults to 'logs/specmap.log'
            max_bytes (int): Maximum size of log file before rotation (default: 10MB)
            backup_count (int): Number of backup files to keep (default: 5)
        """
        self.log_file = log_file
        
        # Create logs directory if it doesn't exist
        log_dir = os.path.dirname(log_file)
        if log_dir and not os.path.exists(log_dir):
            try:
                os.makedirs(log_dir, exist_ok=True)
            except Exception as e:
                # Fallback to current directory if can't create logs directory
                print(f"Warning: Could not create logs directory '{log_dir}': {e}")
                self.log_file = 'specmap.log'
                log_dir = '.'
        
        # Configure logger
        self.logger = logging.getLogger('SpecMap')
        self.logger.setLevel(logging.DEBUG)
        
        # Prevent duplicate handlers if re-initialized
        if not self.logger.handlers:
            # Create rotating file handler
            try:
                file_handler = RotatingFileHandler(
                    self.log_file,
                    maxBytes=max_bytes,
                    backupCount=backup_count,
                    encoding='utf-8'
                )
                file_handler.setLevel(logging.DEBUG)
                
                # Add filter to ensure context attribute exists
                file_handler.addFilter(ContextFilter())
                
                # Create formatter
                formatter = logging.Formatter(
                    '[%(asctime)s] [%(levelname)s] [%(context)s] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S'
                )
                file_handler.setFormatter(formatter)
                
                # Add handler to logger
                self.logger.addHandler(file_handler)
                
                # Also add console handler for development
                console_handler = logging.StreamHandler()
                console_handler.setLevel(logging.INFO)
                console_handler.addFilter(ContextFilter())
                console_handler.setFormatter(formatter)
                self.logger.addHandler(console_handler)
                
            except Exception as e:
                print(f"Error setting up logging: {e}")
                # Create a basic console-only logger as fallback
                console_handler = logging.StreamHandler()
                console_handler.setLevel(logging.INFO)
                console_handler.addFilter(ContextFilter())
                
                # Create formatter
                formatter = logging.Formatter(
                    '[%(asctime)s] [%(levelname)s] [%(context)s] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S'
                )
                console_handler.setFormatter(formatter)
                self.logger.addHandler(console_handler)
        
        # Log initialization
        self.logger.info("ErrorEngine initialized", extra={'context': 'Initialization'})
    
    def error(self, exception, context="Unknown", **kwargs):
        """
        Handle an error by logging details and displaying user-friendly message.
        
        Args:
            exception (Exception): The exception that occurred
            context (str): Description of what was being done when error occurred
            **kwargs: Additional context information (e.g., filename, action, etc.)
        
        Example:
            error_engine.error(
                exception=FileNotFoundError("File not found"),
                context="Loading spectrum data",
                filename="/path/to/file.txt",
                action="read_file"
            )
        """
        # Build context string
        context_parts = [context]
        for key, value in kwargs.items():
            context_parts.append(f"{key}={value}")
        context_str = " | ".join(context_parts)
        
        # Log the full exception with traceback
        exc_type = type(exception).__name__
        exc_msg = str(exception)
        exc_traceback = ''.join(traceback.format_tb(exception.__traceback__))
        
        log_message = f"Exception: {exc_type}: {exc_msg}\nStacktrace:\n{exc_traceback}"
        
        self.logger.error(
            log_message,
            extra={'context': context_str},
            exc_info=False  # We already formatted the traceback
        )
        
        # Display user-friendly error message
        user_message = self._format_user_message(exception, context, kwargs)
        try:
            #messagebox.showerror("Error", user_message)
            print(f"ERROR: {user_message}")
        except Exception as e:
            # If messagebox fails (e.g., no display), just print
            print(f"Error (messagebox failed): {user_message}")
            print(f"Messagebox error: {e}")
    
    def warning(self, message, context="Unknown", **kwargs):
        """
        Handle a warning by logging and optionally displaying to user.
        
        Args:
            message (str): Warning message
            context (str): Description of context
            **kwargs: Additional context information
        
        Example:
            error_engine.warning(
                message="Low signal detected",
                context="Spectrum analysis",
                threshold=100,
                actual=45
            )
        """
        # Build context string
        context_parts = [context]
        for key, value in kwargs.items():
            context_parts.append(f"{key}={value}")
        context_str = " | ".join(context_parts)
        
        # Log the warning
        self.logger.warning(
            message,
            extra={'context': context_str}
        )
        
        # Display warning to user
        user_message = f"{context}\n\n{message}"
        if kwargs:
            details = "\n".join([f"{k}: {v}" for k, v in kwargs.items()])
            user_message += f"\n\nDetails:\n{details}"
        
        try:
            #messagebox.showwarning("Warning", user_message)
            print(f"WARNING: {user_message}")
        except Exception as e:
            # If messagebox fails, just print
            print(f"Warning (messagebox failed): {user_message}")
            print(f"Messagebox error: {e}")
    
    def get_logger(self):
        """
        Get the configured logger instance for direct logging.
        
        Returns:
            logging.Logger: The configured logger
            
        Example:
            logger = error_engine.get_logger()
            logger.info("Processing started", extra={'context': 'Data Processing'})
            logger.debug("Processing file X", extra={'context': 'File Processing'})
        """
        return self.logger
    
    def _format_user_message(self, exception, context, kwargs):
        """
        Format a user-friendly error message from exception details.
        
        Args:
            exception (Exception): The exception
            context (str): Context description
            kwargs (dict): Additional context
            
        Returns:
            str: Formatted user message
        """
        exc_type = type(exception).__name__
        exc_msg = str(exception)
        
        # Create user-friendly message
        message_parts = [
            f"An error occurred: {context}",
            "",
            f"Error type: {exc_type}",
            f"Message: {exc_msg}"
        ]
        
        # Add context details if provided
        if kwargs:
            message_parts.append("")
            message_parts.append("Details:")
            for key, value in kwargs.items():
                message_parts.append(f"  {key}: {value}")
        
        message_parts.append("")
        message_parts.append("Please check the log file for more details:")
        message_parts.append(f"  {self.log_file}")
        
        return "\n".join(message_parts)


# Module-level convenience function for quick initialization
_default_error_engine = None

def get_default_error_engine():
    """
    Get or create the default ErrorEngine instance.
    
    This provides a singleton error engine that can be safely accessed
    from any thread or context in the application. The error engine
    should be initialized once in the main application (main9.py) to
    configure the log file location and rotation settings.
    
    Returns:
        ErrorEngine: The default error engine instance
        
    Note:
        If called before explicit initialization in main9.py, this will
        create a default instance with default settings (logs/specmap.log).
    """
    global _default_error_engine
    if _default_error_engine is None:
        _default_error_engine = ErrorEngine()
    return _default_error_engine


# Example usage demonstration
if __name__ == "__main__":
    # Example 1: Basic usage
    print("Example 1: Basic error handling")
    engine = ErrorEngine()
    
    try:
        # Simulate an error
        raise FileNotFoundError("spectrum_001.txt not found")
    except Exception as e:
        engine.error(
            exception=e,
            context="Reading spectrum file",
            filename="spectrum_001.txt",
            action="file_open"
        )
    
    # Example 2: Warning usage
    print("\nExample 2: Warning handling")
    engine.warning(
        message="Signal intensity below threshold",
        context="Data quality check",
        threshold=100,
        actual_value=45,
        filename="spectrum_002.txt"
    )
    
    # Example 3: Direct logger usage
    print("\nExample 3: Direct logger usage")
    logger = engine.get_logger()
    logger.info("Application initialized successfully", extra={'context': 'Startup'})
    logger.debug("Loading configuration file", extra={'context': 'Configuration'})
    
    print(f"\nLog file created at: {engine.log_file}")
