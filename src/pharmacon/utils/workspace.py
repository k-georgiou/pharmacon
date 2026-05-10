"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

Module :mod:`pharmacon.utils.workspace`.
"""
import os
import atexit
import shutil
import tempfile

from pathlib import Path
from typing import Optional

from pharmacon.logger import get_logger, PharmaconLogger




__all__ = [
    "logger",
    "PharmaconWorkspace",
]


logger: PharmaconLogger = get_logger(__name__)


class PharmaconWorkspace:
    """
    Manages a workspace for Pharmacon, including initialization, working directory,
    temporary directory handling, and cleanup.

    This class provides utilities to manage the current working environment by handling
    a working directory and an optional temporary directory. It ensures proper setup
    and cleanup of resources necessary for workspace operations.

    :ivar origin_cwd: The original current working directory when the object is instantiated.
    :type origin_cwd: Path
    :ivar working_directory: The directory where workspace operations are conducted.
    :type working_directory: Path
    :ivar temp_directory: The optional temporary directory created for workspace operations.
    :type temp_directory: Optional[Path]
    """

    def __init__(self, *, is_tmp_dir_needed: bool, chdir_to_tmp: bool = False,
                 working_directory: str | Path | None = None, temp_directory_base: str | Path | None = None,
                 cleanup_on_exit: bool = True) -> None:
        """
        This class initializes and manages working and temporary directories for a process.
        It supports creating a temporary directory if needed, optionally changing the working
        directory to the temporary one, and performing cleanup actions upon program exit.

        :param is_tmp_dir_needed: Boolean indicating if a temporary directory is required.
        :param chdir_to_tmp: Boolean indicating if the working directory should switch
            to the temporary directory after its creation.
        :param working_directory: Path for the main working directory. If not provided, defaults to
            a value from the PHARMACON_WORKDIR environment variable or the current directory.
        :param temp_directory_base: Base path for the temporary directory. If not specified,
            defaults to a value from the PHARMACON_TEMPDIR environment variable or None.
        :param cleanup_on_exit: Boolean indicating whether to clean up the temporary directory
            upon program exit.
        """

        # Origin
        self.origin_cwd: Path = Path.cwd().resolve()
        logger.debug("Origin CWD: %s", self.origin_cwd)

        # Working directory
        wd = (
            working_directory
            or os.environ.get("PHARMACON_WORKDIR")
            or "."
        )

        self.working_directory: Path = Path(wd).expanduser().resolve()
        self.working_directory.mkdir(parents=True, exist_ok=True)
        logger.debug("Working directory: %s", self.working_directory)

        # Temporary directory (optional)
        self.temp_directory: Optional[Path] = None
        self._cleanup_on_exit: bool = cleanup_on_exit

        if is_tmp_dir_needed:
            base = (
                temp_directory_base
                or os.environ.get("PHARMACON_TEMPDIR")
                or None
            )

            base_path: Path | None = Path(base).expanduser().resolve() if base else None

            if base_path is not None:
                base_path.mkdir(parents=True, exist_ok=True)

            self.temp_directory = Path(
                tempfile.mkdtemp(prefix="pharmacon-", dir=base_path)
            ).resolve()

            logger.info("Temporary directory created: %s", self.temp_directory)

            if cleanup_on_exit:
                atexit.register(self.cleanup)

            if chdir_to_tmp:
                self.chdir_to_temp()

    def chdir_to_temp(self) -> None:
        """
        Changes the current working directory to the temporary directory.

        This method sets the current working directory to the pre-initialized
        temporary directory. If the temporary directory is not initialized,
        an exception will be raised.

        :raises RuntimeError: If the temporary directory is not initialized.
        """
        if not self.temp_directory:
            raise RuntimeError("Temp directory not initialized.")

        if not self.temp_directory.is_dir():
            raise RuntimeError(f"Temp directory does not exist: {self.temp_directory}")

        logger.debug("Changing directory: %s -> %s", Path.cwd(), self.temp_directory)
        os.chdir(self.temp_directory)

    def restore_cwd(self) -> None:
        """
        Restores the current working directory to its original state.

        This method changes the current working directory back to the initial
        directory stored in `self.origin_cwd`. It ensures that any temporary
        changes to the working directory made during the program execution
        are reverted.

        :return: None
        """
        logger.debug("Restoring cwd -> %s", self.origin_cwd)
        os.chdir(self.origin_cwd)

    def cleanup(self) -> None:
        """
        Cleans up the temporary directory if it has been set.

        This function checks if a temporary directory is set and removes
        it if it exists. It also logs the cleanup event. After the cleanup,
        the `temp_directory` attribute is reset to None to indicate that
        the temporary resource has been cleaned up.

        :return: None
        """
        if self.temp_directory and self.temp_directory.is_dir():
            logger.info("Cleaning up temp directory: %s", self.temp_directory)
            shutil.rmtree(self.temp_directory, ignore_errors=True)
            self.temp_directory = None

    def __enter__(self) -> "PharmaconWorkspace":
        """
        Ensures the context management protocol is adhered to when entering
        the runtime context provided by the PharmaconWorkspace instance.

        :return: Returns the current PharmaconWorkspace instance.
        :rtype: PharmaconWorkspace
        """
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        """
        Handles the exit operation for a context manager. This method restores the
        original working directory and performs cleanup if required.

        :param exc_type: The type of exception raised during the context block, or
            None if no exception was raised.
        :param exc_val: The exception instance raised during the context block, or
            None if no exception was raised.
        :param exc_tb: The traceback object associated with the exception raised during
            the context block, or None if no exception was raised.
        :return: None
        """
        self.restore_cwd()
        if self._cleanup_on_exit:
            self.cleanup()

    def __repr__(self) -> str:
        """
        Constructs and returns a string representation of the object for debugging purposes.

        :return: A string representation of the object, including the class name and values of
                 ``working_directory`` and ``temp_directory`` attributes.
        :rtype: str
        """
        return (
            f"{self.__class__.__name__}("
            f"working_directory={self.working_directory!r}, "
            f"temp_directory={self.temp_directory!r})"
        )
