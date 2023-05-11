#!/usr/bin/python3
import numpy.typing as npt
import pandas as pd
import numpy as np

from typing import Union, Dict, List
from collections import defaultdict
from glob import glob
import pathlib
import sys


if sys.version_info < (3, 11):
    StrArray = np.ndarray
else:
    from typing import TypeAlias
    StrArray: TypeAlias = npt.NDArray[np.str_]


class ParseRes:
    def parse(self, file) -> None:
        self.path: pathlib.Path = pathlib.Path(file)
        lines: StrArray = self._get_lines(self.path)
        self.rw: float = self._get_rw(lines)
        lines: StrArray = self._reduce_lines(lines, 15, "--")
        self.res: Dict[str, float] = self._dict_from_lines(lines)
        self.res["rw"] = self.rw
        self.name: str = self.path.stem


    def _get_lines(self, path: pathlib.Path) -> StrArray:
        text: str = path.read_text()
        lines: StrArray = np.array(text.split("\n"))
        return lines


    def _reduce_lines(self, lines: StrArray, lb: int, ubp: str) -> StrArray:
        lines_lb: StrArray = lines[lb:]
        ub: int = [i for i, l in enumerate(lines_lb) if ubp in l][0]
        lines_lub: StrArray = lines_lb[:ub - 2]
        return lines_lub
    

    def _get_rw(self, lines: StrArray) -> float:
        rw_line: str = lines[11]
        rw_line = rw_line.split()[-1]
        rw: float = float(rw_line)
        return rw


    def _dict_from_lines(self, lines: StrArray) -> Dict[str, float]:
        lines_pm = np.char.partition(lines, " +/-")[:, 0]
        name_value_ws = np.char.partition(lines_pm, " ")[:, [0, 2]]
        name_value = np.char.strip(name_value_ws)
        res: Dict[str, float] = {
                name: float(value) for name, value in name_value
                }
        res_dd = defaultdict(lambda: np.nan, res)
        return res_dd

    
class ParseResDir(ParseRes):
    def __init__(self, dir_path: str, par_path: str, filter: str='') -> None:
        self.parse_dir(dir_path, par_path, filter)


    def parse_dir(self, dir_path: str, par_path, filter: str) -> None:
        files: List[str] = glob(dir_path)
        parameter_path: str = par_path
        parameter: StrArray = self.read_parameter(parameter_path)
        if filter != '':
            parameter = np.array([p for p in parameter if not p.startswith(filter)])
        self.df: pd.DataFrame = self.create_df(files, parameter)


    def read_parameter(self, param_path: str) -> StrArray:
        parameter: StrArray = np.loadtxt(param_path, delimiter=",", dtype="str")
        return parameter


    def create_df(
            self,
            files: List[str],
            param: StrArray
            ) -> pd.DataFrame:
    
        data: List[Dict[str, Union[str, float]]] = []
        
        for file in files:
            pr = ParseRes()
            pr.parse(file)
            #if not set(pr.res.keys()).issubset(param):
            #    diff = set(param) - set(pr.res.keys())
            #    raise ValueError(f"fitres contains invalid keys: {diff}")
            _ = [pr.res[p] for p in param]
            comb = pr.res | {'file_name': pr.name}
            comb = defaultdict(lambda: np.nan, sorted(comb.items()))
            data.append(comb)
        
        data_keys: List[str] = list(data[0].keys())
        data_matrix = np.stack([list(d.values()) for d in data], axis=0)
        df = pd.DataFrame(data_matrix, columns=data_keys)
        df = df[param]
        
        for p in param[1:]:
            df[p] = df[p].astype(float)
        
        return df

