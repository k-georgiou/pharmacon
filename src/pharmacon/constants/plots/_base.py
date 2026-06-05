"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Shared foundation for the ``pharmacon.constants.plots`` package:

- global validation constants (font lists, line styles, legend locations),
- the ``PlotSettingsBase`` dataclass-mixin that every plot settings class
  inherits from (alias registry, ``from_dict``, ``_safe_*`` coercion
  helpers, the ``_warn`` bookkeeper).

Concrete plot classes live in sibling modules (``pli``, ``ppi``,
``universal``, ``pca``) and are all re-exported from the package
``__init__``.
"""
import re
import warnings
from argparse import Namespace
from matplotlib import font_manager
from dataclasses import fields, field
from dataclasses import dataclass, asdict
from matplotlib.colors import is_color_like
from typing import ClassVar, Dict, Type, Tuple, Any, Final, Set, List

from pharmacon.logger import get_logger


__all__ = [
    "VALID_EXTENSIONS",
    "VALID_FONT_WEIGHTS",
    "AVAILABLE_FONTS",
    "VALID_LINE_STYLES",
    "VALID_LEGEND_LOCS",
    "PlotSettingsBase",
    "logger",
    # re-exported helpers the subclasses rely on
    "Namespace",
    "dataclass",
    "field",
    "fields",
    "asdict",
    "ClassVar",
    "Tuple",
    "List",
    "Final",
    "Set",
]


warnings.filterwarnings("ignore")


logger = get_logger(__name__)


VALID_EXTENSIONS: Final[Set] = {".png", ".jpg", ".jpeg", ".svg", ".pdf", ".tif", ".tiff"}

VALID_FONT_WEIGHTS: Final[Set] = {"ultralight", "light", "normal", "regular", "book",
                                  "medium", "roman", "semibold", "demibold", "demi",
                                  "bold", "heavy", "extra bold", "black"}
AVAILABLE_FONTS: Final[Set] = {f.name.lower() for f in font_manager.fontManager.ttflist}

VALID_LINE_STYLES: Final[Set] = {"-", "--", "-.", ":", "solid", "dashed", "dashdot", "dotted"}

VALID_LEGEND_LOCS: Final[Set] = {"best", "upper right", "upper left", "lower left", "lower right",
                                 "right", "center left", "center right", "lower center", "upper center",
                                 "center"}


class PlotSettingsBase:
    """
    Base class for managing plot settings with validation and safe coercion.

    This class serves as a foundation for subclasses that handle plot settings.
    It includes a validation mechanism, safe data type coercion methods, and a
    registry system for managing class aliases. Subclasses should implement the
    `_validate_fields` method to define specific validation logic.

    :ivar alias: Tuple of aliases for this class. Used for registry purposes.
    :type alias: Tuple[str, ...]
    """

    _registry: ClassVar[Dict[str, Type["PlotSettingsBase"]]] = {}
    alias: ClassVar[Tuple[str, ...]] = ()

    # CLASS REGISTRATION
    def __init_subclass__(cls, **kwargs):
        """
        Initialize subclass with optional alias registration.

        This method ensures that subclasses can register themselves with
        aliases in the class registry. It raises an exception if duplicate
        aliases are detected during registration.
        """
        super().__init_subclass__(**kwargs)

        for name in getattr(cls, "alias", ()):
            key = name.upper()
            if key in cls._registry:
                raise ValueError(f"Duplicate alias detected: {key}")
            cls._registry[key] = cls

    def _warn(self, message: str):
        """
        Logs a warning message and increments the internal warnings counter.

        :param message: The warning message to log.
        :type message: str
        """
        logger.warning(message)
        self._current_warnings += 1

    # CONSTRUCTION
    @classmethod
    def from_dict(cls, overrides: dict):
        """
        Creates an instance of the class using a dictionary of overrides. The dictionary
        should contain keys corresponding to the field names in the class, and its values
        will be used to set the class fields if they match valid fields. After assigning
        fields, the instance validates itself by invoking the class's ``validate`` method.

        :param overrides: A dictionary containing keys and values where the keys match
            valid class fields. These fields will be overridden with the provided values.
        :type overrides: dict
        :return: A new instance of the class with the overridden values applied.
        :rtype: cls
        """
        instance = cls()

        valid_fields = {f.name for f in fields(cls)}

        for key, value in overrides.items():
            if key in valid_fields:
                setattr(instance, key, value)

        instance.validate()
        return instance

    # VALIDATION ENTRYPOINT
    def validate(self) -> None:
        """
        Validates fields in the object and resets the warning count.

        This method performs validation on the fields associated with the object,
        ensuring that they meet predefined criteria. Before running the validation,
        it resets the count of current warnings.

        :return: None
        """
        self._current_warnings = 0
        self._validate_fields()

    def _validate_fields(self) -> None:
        """
        Subclasses override this.
        Must:
            - normalize values
            - assign defaults
            - call self._warn(...) instead of raising
        """
        pass

    # SAFE COERCION HELPERS
    def _safe_int(self, value: Any, default: int,
                  min_val=None, max_val=None) -> int:
        """
        Converts a value to an integer safely while handling edge cases and ensuring it
        falls within the specified range limits if provided.

        This method attempts to convert a value to an integer, and if unsuccessful,
        it will return a default value. Additionally, the method validates whether the
        converted value falls within the bounds specified by `min_val` and `max_val`.
        If the value is out of bounds, the default value is returned.

        :param value: The value to be converted to an integer. Can be of any type.
        :param default: The fallback integer to use if conversion fails or the value
            does not satisfy bounds.
        :param min_val: Optional. Minimum bound for the converted integer. If None,
            no lower bound is enforced.
        :param max_val: Optional. Maximum bound for the converted integer. If None,
            no upper bound is enforced.
        :return: The resulting integer value if conversion succeeds and bounds are
            satisfied, or the `default` value otherwise.
        :rtype: int
        """
        try:
            value = int(value)
        except Exception:
            self._warn(f"Invalid int '{value}', using default {default}")
            return default

        if min_val is not None and value < min_val:
            self._warn(f"Int '{value}' < {min_val}, using default {default}")
            return default

        if max_val is not None and value > max_val:
            self._warn(f"Int '{value}' > {max_val}, using default {default}")
            return default

        return value

    def _safe_float(self, value: Any, default: float,
                    min_val=None, max_val=None) -> float:
        """
        Converts the input value to a float while ensuring it falls within the
        specified range. If the conversion fails or the value falls outside the
        range, a default value is returned instead. Logs warnings for invalid
        values or range violations.

        :param value: The value to be converted to a float.
        :param default: The default float value to return in case of an invalid
            conversion or if the value is out of range.
        :param min_val: The minimum allowable value for the converted float. If
            the value is less than `min_val`, the default is returned. Defaults to None.
        :param max_val: The maximum allowable value for the converted float. If
            the value is greater than `max_val`, the default is returned. Defaults to None.
        :return: The float value if conversion is successful and within the
            specified range; otherwise, the default value.
        :rtype: float
        """
        try:
            value = float(value)
        except Exception:
            self._warn(f"Invalid float '{value}', using default {default}")
            return default

        if min_val is not None and value < min_val:
            self._warn(f"Float '{value}' < {min_val}, using default {default}")
            return default

        if max_val is not None and value > max_val:
            self._warn(f"Float '{value}' > {max_val}, using default {default}")
            return default

        return value

    def _safe_bool(self, value: Any, default: bool) -> bool:
        """
        Converts a given value into a boolean, providing a fallback to a default
        value if the conversion is not possible. Validations include conversions
        from boolean, string representations, and ensuring logical defaults
        for invalid inputs.

        :param value: The input value to evaluate and convert to a boolean.
        :type value: Any
        :param default: The default boolean value to fall back on if the given
            value cannot be interpreted as a boolean.
        :return: The resultant boolean value based on conversion or default.
        :rtype: bool
        """
        if isinstance(value, bool):
            return value

        # Accept native int 0 / 1 (e.g. "flag = 0" from an INI file).
        if isinstance(value, int) and value in (0, 1):
            return bool(value)

        if isinstance(value, str):
            v = value.strip().lower()
            if v in {"true", "1", "yes", "y"}:
                return True
            if v in {"false", "0", "no", "n"}:
                return False

        self._warn(f"Invalid bool '{value}', using default {default}")
        return default

    def _safe_color(self, value: Any, default: str) -> str:
        """
        Ensures that the provided value is a valid color. If it is not valid, a default
        color is returned instead, and a warning is issued. This function is intended
        to handle cases where color validation is necessary.

        :param value: The input value to validate as a color.
        :type value: Any
        :param default: The default color value to use if the input value is not a valid color.
        :type default: str
        :return: The validated color if valid, otherwise the default color.
        :rtype: str
        """
        value = str(value).strip()
        if is_color_like(value):
            return value

        self._warn(f"Invalid color '{value}', using default '{default}'")
        return default

    # REGISTRY HELPERS
    @classmethod
    def resolve(cls, alias: str) -> Type["PlotSettingsBase"]:
        """
        Resolves an alias to the corresponding class type from the registry.

        This class method is used to look up and return a class type corresponding to the
        provided alias. The alias is matched in a case-insensitive manner by converting it
        to uppercase before lookup. It ensures dynamic resolution of class types registered
        under specific aliases.

        :param alias: A string identifier used to reference a specific class type. The
            alias is matched in a case-insensitive manner.
        :type alias: str
        :return: The class type associated with the provided alias.
        :rtype: Type["PlotSettingsBase"]
        :raises KeyError: If the alias does not exist in the registry.
        """
        return cls._registry[alias.upper()]

    @classmethod
    def get_registry(cls):
        """
        Retrieve the current state of the registry as a dictionary.

        This class method provides access to the registry, which is
        maintained as a mapping in the class-level storage. It returns
        a copy of the registry to ensure that direct modifications cannot
        occur to the underlying data.

        :return: A dictionary containing the registered items.
        :rtype: dict
        """
        return dict(cls._registry)

    @classmethod
    def get_all_aliases(cls):
        """
        Provides a class-level utility method to retrieve all registered aliases
        associated with the `cls` class. Aliases are used to map and identify
        specific entries in the `_registry` dictionary of the class.

        :raises KeyError: If the `_registry` attribute is not defined in the
            class.
        :rtype: list
        :return: A list of alias names registered within the class.
        """
        return list(cls._registry.keys())


