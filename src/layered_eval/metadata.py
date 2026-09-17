"""Strict cell metadata used by both counts preparation and bundle validation."""
import numpy as np
import pandas as pd


def validate_cell_metadata(cells):
    """Return normalized metadata without treating strings or missing flags as truthy.

    CSV boolean fields accept only true/false (case-insensitive, outer whitespace
    ignored). Numeric 0/1, empty values and other spellings are rejected. Both
    flags and the original source count must be constant within a guide.
    """
    required = {'cell_id', 'guide_id', 'target_id', 'is_control',
                'source_excluded_guide', 'source_guide_n'}
    missing = sorted(required.difference(cells.columns))
    if missing:
        raise ValueError('cells.csv missing required metadata column(s): ' + ', '.join(missing))
    cells = cells.copy()
    for column in ('cell_id', 'guide_id', 'target_id'):
        if cells[column].isna().any() or cells[column].astype(str).str.strip().eq('').any():
            raise ValueError(f'cells.csv {column} must be nonmissing and nonempty')
    if not cells.cell_id.is_unique:
        raise ValueError('cell IDs must be unique and nonmissing')
    for column in ('is_control', 'source_excluded_guide'):
        normalized = []
        for position, value in enumerate(cells[column]):
            if isinstance(value, (bool, np.bool_)):
                normalized.append(bool(value))
            elif isinstance(value, str) and value.strip().lower() in ('true', 'false'):
                normalized.append(value.strip().lower() == 'true')
            else:
                raise ValueError(
                    f'cells.csv {column} requires explicit true/false at data row '
                    f'{position + 1} (cell_id={cells.cell_id.iloc[position]!r}); '
                    'missing values, 0/1 and other spellings are not accepted')
        cells[column] = np.asarray(normalized, dtype=bool)
    count = pd.to_numeric(cells.source_guide_n, errors='coerce').to_numpy(dtype=float)
    if not np.isfinite(count).all() or (count < 1).any() or (count != np.floor(count)).any():
        raise ValueError('cells.csv source_guide_n must contain positive integer original guide counts')
    cells['source_guide_n'] = count.astype(np.int64)
    for column in ('is_control', 'source_excluded_guide', 'source_guide_n'):
        inconsistent = cells.groupby('guide_id', sort=False)[column].nunique()
        bad = inconsistent[inconsistent > 1].index.tolist()
        if bad:
            raise ValueError(
                f'cells.csv {column} is inconsistent within guide_id {bad[:5]!r}; '
                'resolve the source annotation before preparing or validating a bundle')
    return cells
