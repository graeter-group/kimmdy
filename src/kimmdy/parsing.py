from pathlib import Path
from collections.abc import Iterable
from typing import Generator

def is_comment(l : str):
    return(len(l) == 0 or l[0] in ['#', '\n'])

def get_sections(seq : Iterable[str], section_marker : str) -> Generator[list[str], None, None]:
    data = []
    for line in seq:
        if is_comment(line): continue
        if line.startswith(section_marker):
            if data:
                yield data
                data = []
        data.append(line.strip('\n'))
    if data:
        yield data

def read_topol(path : Path) :
    with open(path, 'r') as f:
        sections = get_sections(f, '[')
        return { title.strip('[] \n'): content for title, *content in sections }

def write_topol(d : dict[str, str], outfile : Path):
    with open(outfile, 'w') as f:
        for key, value in d.items():
            s = f'[ {key} ]\n'
            f.write(s)
            s = '\n'.join(value) + '\n\n'
            f.write(s)

