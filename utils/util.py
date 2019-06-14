# util.py

from pathlib import Path
from typing import Any, Dict, NoReturn, Tuple, Union


def subset_dict(key_range: Tuple[int, int], original_dict: Dict[int, Any]) -> Dict[int, Any]:
    """
    Find the subset of keys within the given key range. Key must be integer values.

    Parameters
    ----------
    key_range : tuple of int
        A tuple pair representing the boundaries of key values within the subset.

    original_dict : dict of {int, Any}
        The full dictionary of integer keys to any value.

    Returns
    -------
    dict of {int, Any}
        A subset dictionary with only entries where the keys lie within the key range bounds.
    """
    subset_dict: Dict[int, Any] = dict(
        (k, original_dict[k]) for k in original_dict.keys() if key_range[0] <= k < key_range[1]
    )

    return subset_dict


def reverse_tuple_pair(tup: Tuple[Any, Any]) -> Union[NoReturn, Tuple[Any, Any]]:
    """
    Reverse a tuple pair. Pairs only.

    Parameters
    ----------
    tup : tuple of any values
        A tuple pair of any values.

    Returns
    -------
    tuple
        The reversed tuple pair.

    Raises
    ------
    ValueError
        If the tuple isn't a pair.
    """
    if len(tup) != 2:
        raise ValueError("Provided tuple isn't a pair.")

    return (tup[1], tup[0])


def ensure_dir(dirpath: Path) -> None:
    """
    Checks if a directory exists, creates it if it doesn't.

    Parameters
    ----------
    dirpath : pathlib.Path
        Path of the directory in question.
    """
    if not dirpath.is_dir():
        dirpath.mkdir(parents=True)
