"""Pharmacon — Molecular Dynamics Suite, developed by Kyriakos Georgiou, 2026.

Shared validation helpers for CLI subcommands.

All functions raise :class:`~pharmacon.command_line.exceptions.ValidationError`
on invalid input.
"""

from __future__ import annotations

from pathlib import Path

from pharmacon.command_line.exceptions import ValidationError


# Path helpers




__all__ = [
    "normalize_path",
    "validate_existing_input_file",
    "validate_output_file",
    "normalize_selection",
    "validate_string_list",
    "validate_non_negative_int",
    "validate_positive_int",
    "validate_bool_flag",
    "validate_logging_level",
    "validate_frame_range",
    "validate_reference_frame",
]


def normalize_path(value: object, field_name: str) -> Path:
    """
    Normalizes a given path-like value into a `Path` object, ensuring validation.
    This function checks that the input is either a string or a `Path` object,
    strips any leading or trailing whitespace, and verifies that it does not
    contain control characters. If the input passes all validation rules, it
    returns the expanded user path as a `Path` object.

    :param value: The input path-like value to normalize.
    :type value: object
    :param field_name: The name of the field used in error messages for
        reference.
    :type field_name: str
    :raises ValidationError: If the value is not path-like, is empty, or
        contains invalid characters.
    :return: The normalized path as a `Path` object.
    :rtype: Path
    """
    if not isinstance(value, (str, Path)):
        raise ValidationError(
            f"Argument '{field_name}' must be a path-like value, "
            f"got {type(value).__name__}."
        )
    text: str = str(value).strip()
    if not text:
        raise ValidationError(f"Argument '{field_name}' cannot be empty.")
    if any(char in text for char in ("\x00", "\r", "\n")):
        raise ValidationError(
            f"Argument '{field_name}' contains invalid control characters."
        )
    return Path(text).expanduser()


def validate_existing_input_file(value: object,
                                 field_name: str,
                                 supported_formats: tuple[str, ...] | list[str] | set[str] | frozenset[str]) -> Path:
    """
    Validates that the provided value corresponds to an existing, readable input file
    with a format explicitly listed in the provided supported formats. It ensures
    the path exists, is a regular file, matches one of the supported formats, and
    can be opened for reading. If any of these conditions fail, a ValidationError
    is raised.

    :param value: The input value expected to be a valid file path.
    :type value: object
    :param field_name: The name of the field related to the file, used for error
        message customization.
    :param supported_formats: A collection of supported file extensions, such as
        (".txt", ".csv"). Formats are case-insensitive and may optionally start
        with a leading dot (e.g., ".ext" or "ext").
    :type supported_formats: tuple[str, ...] | list[str] | set[str] | frozenset[str]
    :raises ValidationError: If the file does not exist, is not a regular file,
        has an unsupported format, or cannot be read.
    :return: The absolute, resolved path to the validated input file.
    :rtype: Path
    """
    path: Path = normalize_path(value, field_name)
    if not path.exists():
        raise ValidationError(
            f"{field_name.capitalize()} file does not exist: '{path}'."
        )
    if not path.is_file():
        raise ValidationError(
            f"{field_name.capitalize()} path is not a regular file: '{path}'."
        )
    supported_suffixes: set[str] = {
        fmt.lower() if str(fmt).startswith(".") else f".{str(fmt).lower()}"
        for fmt in supported_formats
    }
    suffix: str = path.suffix.lower()
    if suffix not in supported_suffixes:
        supported_text: str = ", ".join(sorted(supported_suffixes))
        raise ValidationError(
            f"Unsupported {field_name} format for '{path.name}'. "
            f"Supported formats: {supported_text}."
        )
    try:
        with path.open("rb"):
            pass
    except OSError as exc:
        raise ValidationError(
            f"{field_name.capitalize()} file is not readable: '{path}'."
        ) from exc
    return path.resolve()


def validate_output_file(value: object,
                         field_name: str,
                         *,
                         overwrite: bool,
                         allowed_suffixes: tuple[str, ...] | list[str] | set[str] | None = None
                         ) -> Path:
    """
    Validates the specified output file path for various conditions. It checks if the path is a valid and writable
    file path. The validation includes checks for file existence, allowed suffixes, directory existence,
    and whether overwriting is permitted. Also ensures that the file is not a directory and its parent path
    is a valid directory.

    :param value: The input value representing the file path to validate.
    :type value: object
    :param field_name: The name of the parameter being validated.
    :type field_name: str
    :param overwrite: Determines if the existing file, if any, can be overwritten.
    :type overwrite: bool
    :param allowed_suffixes: A collection of allowed suffixes for the file. These suffixes are normalized to
        lowercase and prefixed with a dot (e.g., '.txt'). If set to None, this validation is skipped.
    :type allowed_suffixes: tuple[str, ...] | list[str] | set[str] | None
    :return: The resolved and validated file path.
    :rtype: Path
    """
    path: Path = normalize_path(value, field_name)
    if path.name in {"", ".", ".."}:
        raise ValidationError(f"Argument '{field_name}' has an invalid file name.")
    if allowed_suffixes is not None:
        normalized_suffixes: set[str] = {
            fmt.lower() if str(fmt).startswith(".") else f".{str(fmt).lower()}"
            for fmt in allowed_suffixes
        }
        suffix: str = path.suffix.lower()
        if suffix not in normalized_suffixes:
            supported_text: str = ", ".join(sorted(normalized_suffixes))
            raise ValidationError(
                f"Unsupported {field_name} format for '{path.name}'. "
                f"Supported formats: {supported_text}."
            )
    parent: Path = path.parent.expanduser()
    if not parent.exists():
        raise ValidationError(
            f"Parent directory for '{field_name}' does not exist: '{parent}'."
        )
    if not parent.is_dir():
        raise ValidationError(
            f"Parent path for '{field_name}' is not a directory: '{parent}'."
        )
    if path.exists() and path.is_dir():
        raise ValidationError(
            f"Argument '{field_name}' points to a directory, not a file: '{path}'."
        )
    if path.exists() and not overwrite:
        raise ValidationError(
            f"{field_name.capitalize()} file already exists: '{path}'. "
            "Use '--overwrite' to replace it."
        )
    try:
        if path.exists():
            with path.open("ab"):
                pass
        else:
            with path.open("xb"):
                pass
            path.unlink()
    except FileExistsError as exc:
        raise ValidationError(
            f"{field_name.capitalize()} file already exists: '{path}'. "
            "Use '--overwrite' to replace it."
        ) from exc
    except OSError as exc:
        raise ValidationError(
            f"Cannot write to '{field_name}' path: '{path}'."
        ) from exc
    return path.resolve()


# Selection helpers
def normalize_selection(value: object,
                        field_name: str,
                        *,
                        required: bool,
                        default: str | None = None
                        ) -> str | None:
    """
    Normalizes the input string and validates it according to the given rules.
    If the value is None, it checks whether the field is required and uses the
    default value if provided. It ensures the string is non-blank and does not
    contain prohibited control characters such as null, carriage return, or
    newline.

    :param value: The input object to be normalized and validated.
    :param field_name: The name of the field to include in error messages.
    :param required: Indicates whether the field is mandatory. If True, a non-empty
        string must be provided unless a default value is specified.
    :param default: The default value to use when the input is None or blank,
        applicable only when the field is not required or a non-blank default is
        provided.
    :return: The normalized string if valid, or the default value if applicable,
        otherwise None.
    :raises ValidationError: If the field is required and the value is invalid or
        missing, or if the value contains invalid characters.
    """
    if value is None:
        if required and default is None:
            raise ValidationError(
                f"Argument '{field_name}' requires a non-empty selection string."
            )
        return default
    if not isinstance(value, str):
        raise ValidationError(
            f"Argument '{field_name}' must be a string, "
            f"got {type(value).__name__}."
        )
    normalized: str = value.strip()
    if not normalized:
        if required and default is None:
            raise ValidationError(f"Argument '{field_name}' cannot be blank.")
        return default
    if any(char in normalized for char in ("\x00", "\r", "\n")):
        raise ValidationError(
            f"Argument '{field_name}' contains invalid control characters."
        )
    return normalized


def validate_string_list(value: object, field_name: str) -> list[str]:
    """
    Validates that the provided value is a non-empty list of non-empty strings. Each string
    in the list is trimmed of surrounding whitespace before being added to the result. If
    the provided value is not a valid list of strings or is empty, or contains any invalid
    strings, a ValidationError is raised.

    :param value: The value to be validated.
    :type value: object
    :param field_name: The name of the field being validated to identify validation errors.
    :type field_name: str
    :return: A list of non-empty strings, each trimmed of surrounding whitespace.
    :rtype: list[str]
    :raises ValidationError: If the value is not a list, is empty, contains invalid
        strings, or strings that are empty or consist solely of whitespace.
    """
    if not isinstance(value, list) or not value:
        raise ValidationError(
            f"Argument '{field_name}' requires at least one value."
        )
    result: list[str] = []
    for i, item in enumerate(value):
        if not isinstance(item, str) or not item.strip():
            raise ValidationError(
                f"{field_name} at position {i + 1} is empty or invalid."
            )
        result.append(item.strip())
    return result


# Scalar validators
def validate_non_negative_int(value: object, field_name: str) -> int:
    """
    Validates that the given value is a non-negative integer. This function ensures
    the value is of type integer, is not a boolean (as booleans are treated as
    integers in Python), and is greater than or equal to 0. If any of these
    conditions are not met, a ValidationError is raised.

    :param value: The value to be validated.
    :type value: object
    :param field_name: The name of the field associated with the value, used
        in the error message for better clarity.
    :type field_name: str
    :return: The validated integer value if it passes all checks.
    :rtype: int
    :raises ValidationError: If the value is not an integer, is a boolean,
        or is a negative integer.
    """
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValidationError(f"Argument '{field_name}' must be an integer.")
    if value < 0:
        raise ValidationError(
            f"Argument '{field_name}' must be >= 0, got {value}."
        )
    return value


def validate_positive_int(value: object, field_name: str) -> int:
    """
    Validates that the provided value is a positive integer. This includes checking that
    the value is of the correct type (int) and greater than zero. If validation fails,
    a `ValidationError` is raised with an appropriate error message.

    :param value: The value to validate. Must be an integer and greater than zero.
    :type value: object
    :param field_name: The name of the field being validated, included in error
        messages for better context.
    :type field_name: str
    :return: The validated positive integer value.
    :rtype: int
    :raises ValidationError: If `value` is not an integer, is a boolean value
        (since booleans are a subclass of integers), or if `value` is not greater
        than zero.
    """
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValidationError(f"Argument '{field_name}' must be an integer.")
    if value <= 0:
        raise ValidationError(
            f"Argument '{field_name}' must be > 0, got {value}."
        )
    return value


def validate_bool_flag(value: object, field_name: str) -> bool:
    """
    Validates that a given value is of type `bool`. This function is useful for ensuring that an argument
    provided to a method or function is a boolean flag. If the value is not a boolean, a `ValidationError`
    is raised to enforce the requirement.

    :param value: The value to be validated as a boolean flag.
    :type value: object
    :param field_name: The name of the field or parameter being validated.
    :type field_name: str
    :return: Returns the validated boolean value if the validation is successful.
    :rtype: bool
    :raises ValidationError: If the provided value is not of type `bool`.
    """
    if not isinstance(value, bool):
        raise ValidationError(f"Argument '{field_name}' must be a boolean flag.")
    return value


def validate_logging_level(value: object, field_name: str) -> str:
    """
    Validates that the given logging level value is a valid string and matches one of the
    accepted logging levels. If the value is not valid, raises a ValidationError.

    :param value: The logging level string to validate.
    :type value: object
    :param field_name: The name of the field corresponding to the logging level value,
        used for error message context.
    :type field_name: str
    :return: The normalized and validated logging level string.
    :rtype: str
    :raises ValidationError: If the value is not a string or not one of the valid logging
        levels.
    """
    if not isinstance(value, str):
        raise ValidationError(
            f"Argument '{field_name}' must be a logging level string."
        )
    normalized: str = value.strip().upper()
    valid_levels: set[str] = {
        "TRACE", "DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL",
    }
    if normalized not in valid_levels:
        valid_text: str = ", ".join(sorted(valid_levels))
        raise ValidationError(
            f"Invalid logging level for '{field_name}': '{value}'. "
            f"Valid levels: {valid_text}."
        )
    return normalized


# Composite validators
def validate_frame_range(data: dict[str, object]) -> tuple[int, int | None, int]:
    """
    Validates and processes the frame range parameters within the provided dictionary.
    Ensures the parameters 'begin', 'end', and 'step' adhere to specified constraints such as
    non-negativity and logical relationships between 'begin' and 'end'.

    :param data: A dictionary containing the keys 'begin', 'end', and 'step', which
        respectively represent the start of the range, optional end of the range, and
        increment for the range.
    :type data: dict[str, object]

    :return: A tuple containing three elements:
        - An integer representing the validated start of the range ('begin').
        - Either an integer or None representing the validated end of the range ('end').
        - An integer representing the validated step increment ('step').
    :rtype: tuple[int, int | None, int]

    :raises ValidationError: If 'end' is less than 'begin' after validation or if any of
        the numerical validations fail for 'begin', 'end', or 'step'.
    """
    begin = validate_non_negative_int(data.get("begin", 0), "begin")
    raw_end: object = data.get("end")
    if raw_end is None:
        end: int | None = None
    else:
        end = validate_non_negative_int(raw_end, "end")
    step = validate_positive_int(data.get("step", 1), "step")
    if end is not None and end < begin:
        raise ValidationError(
            f"Argument 'end' ({end}) must be greater than or equal to "
            f"'begin' ({begin})."
        )
    data["begin"] = begin
    data["end"] = end
    data["step"] = step
    return begin, end, step


def validate_reference_frame(data: dict[str, object],
                             begin: int,
                             end: int | None,
                             field_name: str = "reference_frame"
                             ) -> int:

    """
    Validates the reference frame from the provided data dictionary and ensures it
    falls within the specified range defined by `begin` and `end`. Updates the
    `field_name` key in the `data` dictionary with the validated reference frame and
    returns it.

    :param data: Dictionary containing data to validate.
    :param begin: The minimum valid value for the reference frame.
    :param end: The maximum valid value for the reference frame. If `None`, there is
        no upper limit.
    :param field_name: The key name in the `data` dictionary representing the
        reference frame. Default is "reference_frame".
    :return: The validated non-negative reference frame value.
    :raises ValidationError: If the reference frame is less than `begin` or greater
        than `end` (if `end` is not `None`).
    """
    ref = validate_non_negative_int(data.get(field_name, 0), field_name)
    if ref < begin:
        raise ValidationError(
            f"Reference frame ({ref}) must be >= begin ({begin})."
        )
    if end is not None and ref > end:
        raise ValidationError(
            f"Reference frame ({ref}) must be <= end ({end})."
        )
    data[field_name] = ref
    return ref
