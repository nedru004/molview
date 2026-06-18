"""Convert py3dmol-style selection specs to Mol* interactivity configs."""

from __future__ import annotations

import re
from typing import Any


def _normalize_chain(chain: str | list) -> str | list[str]:
    if isinstance(chain, (list, tuple)):
        return [str(c) for c in chain]
    return str(chain)


def _parse_resi_values(resi: int | str | list) -> dict[str, Any]:
    """Parse py3dmol resi field into Mol* schema fragments."""
    if isinstance(resi, int):
        return {'auth_seq_id': resi}

    if isinstance(resi, str):
        resi = [resi]

    if not isinstance(resi, (list, tuple)):
        return {}

    ids: list[int] = []
    ranges: list[tuple[int, int]] = []

    for item in resi:
        if isinstance(item, int):
            ids.append(item)
            continue

        text = str(item).strip()
        if not text:
            continue

        range_match = re.match(r'^(-?\d+)\s*[-:]\s*(-?\d+)$', text)
        if range_match:
            start, end = int(range_match.group(1)), int(range_match.group(2))
            if start <= end:
                ranges.append((start, end))
            else:
                ranges.append((end, start))
            continue

        if text.isdigit() or (text.startswith('-') and text[1:].isdigit()):
            ids.append(int(text))

    result: dict[str, Any] = {}
    if ids:
        result['auth_seq_id'] = ids if len(ids) > 1 else ids[0]
    if ranges:
        result['ranges'] = ranges
    return result


def _normalize_resn(resn: str | list) -> str | list[str]:
    if isinstance(resn, (list, tuple)):
        return [str(r).upper() for r in resn]
    return str(resn).upper()


def selection_to_elements(sel: dict | None) -> dict | None:
    """
    Convert a py3dmol AtomSelectionSpec dict to a Mol* StructureElement.Schema.

    Supported fields: chain, resi, resn, elem, atom, serial, hetflag, model.
    """
    if not sel:
        return None

    elements: dict[str, Any] = {}
    items: dict[str, Any] = {}

    if 'chain' in sel:
        chain = _normalize_chain(sel['chain'])
        if isinstance(chain, list):
            items['auth_asym_id'] = chain
        else:
            elements['auth_asym_id'] = chain

    if 'resi' in sel:
        resi_data = _parse_resi_values(sel['resi'])
        if 'auth_seq_id' in resi_data:
            value = resi_data['auth_seq_id']
            if isinstance(value, list):
                items['auth_seq_id'] = value
            else:
                elements['auth_seq_id'] = value
        if 'ranges' in resi_data:
            elements['ranges'] = resi_data['ranges']

    if 'resn' in sel:
        resn = _normalize_resn(sel['resn'])
        if isinstance(resn, list):
            items['auth_comp_id'] = resn
        else:
            elements['auth_comp_id'] = resn

    if 'elem' in sel:
        elem = sel['elem']
        if isinstance(elem, (list, tuple)):
            items['type_symbol'] = [str(e).upper() for e in elem]
        else:
            elements['type_symbol'] = str(elem).upper()

    if 'atom' in sel:
        atom = sel['atom']
        if isinstance(atom, (list, tuple)):
            items['auth_atom_id'] = [str(a) for a in atom]
        else:
            elements['auth_atom_id'] = str(atom)

    if 'serial' in sel:
        serial = sel['serial']
        if isinstance(serial, (list, tuple)):
            items['id'] = [int(s) for s in serial]
        else:
            elements['id'] = int(serial)

    if 'hetflag' in sel:
        elements['group_PDB'] = 'HETATM' if sel['hetflag'] else 'ATOM'

    if 'model' in sel:
        elements['model_num'] = int(sel['model'])

    if items:
        elements['items'] = items

    return elements if elements else None


def normalize_style(style: dict | None) -> dict:
    """Extract color and representation hints from a py3dmol style dict."""
    if not style:
        return {}

    result: dict[str, Any] = {}

    for repr_type in ('cartoon', 'stick', 'sphere', 'line'):
        if repr_type not in style:
            continue

        repr_style = style[repr_type]
        if not isinstance(repr_style, dict):
            result[repr_type] = True
            continue

        entry: dict[str, Any] = {}
        color = repr_style.get('color')
        if color:
            entry['color'] = _normalize_color(color)
        if repr_style:
            entry['params'] = repr_style
        result[repr_type] = entry if entry else True

    return result


def _normalize_color(color: str) -> str:
    color = str(color).strip()
    if color.startswith('0x') or color.startswith('0X'):
        return '#' + color[2:].zfill(6)
    if not color.startswith('#'):
        named = {
            'white': '#FFFFFF',
            'black': '#000000',
            'red': '#FF0000',
            'green': '#00FF00',
            'blue': '#0000FF',
            'yellow': '#FFFF00',
            'purple': '#800080',
            'orange': '#FFA500',
            'cyan': '#00FFFF',
            'magenta': '#FF00FF',
            'gray': '#808080',
            'grey': '#808080',
        }
        return named.get(color.lower(), f'#{color}')
    return color


def build_interaction_spec(
    sel: dict | None,
    action: str,
    *,
    color: str | None = None,
    style: dict | None = None,
    additive: bool = False,
) -> dict:
    """Build a JSON-serializable interaction spec for the viewer template."""
    spec: dict[str, Any] = {
        'action': action,
        'additive': additive,
    }

    elements = selection_to_elements(sel)
    if elements:
        spec['elements'] = elements

    if color:
        spec['color'] = _normalize_color(color)

    if style:
        spec['style'] = normalize_style(style)

    if sel and sel.get('invert'):
        spec['invert'] = True

    return spec
