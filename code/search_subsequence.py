import re
import subprocess
import numpy as np

from time import time
from sortedcontainers import SortedDict

OUPUT_BASE_REGEXP = \
    r"\r\nLocation\s+:\s+(?P<location>\w+)\r\n" + \
    r"Distance\s+:\s+(\d+\.\d+)\r\n" + \
    r"Data Scanned\s+:\s+(\w+)\r\n" + \
    r"Total Execution Time\s+:\s+(\d+\.\d+) sec\r\n\r\n" + \
    r"Pruned by LB_Kim\s+:\s+(\d+\.\d+)%\r\n" + \
    r"Pruned by LB_Keogh\s+:\s+(\d+\.\d+)%\r\n" + \
    r"Pruned by LB_Keogh2\s+:\s+(\d+\.\d+)%\r\n" + \
    r"DTW Calculation\s+:\s+(\d+\.\d+)%\r\n" + \
    r"\r\nlocate\s+dist\r\n----------------\r\n"

OUPUT_STARTS_REGEXP = r"(\d+)\s+(\d+\.\d+)\r\n"

BIN_EXE = r"bin/UCR_MDTW.exe"
BIN_EXE_NO_OPTIM = r"bin/UCR_MDTW_no_optim.exe"
MAX_DIST = 1e20


def zscore(a, axis, ddof):
    a = np.asanyarray(a)
    mns = a.mean(axis=axis)
    sstd = a.std(axis=axis, ddof=ddof)

    # very bad lifehack
    sstd[sstd == 0] = 1

    return (a - mns) / sstd


def search_dtw(closest_series_num: int, subseq_len: int, warp_window: float, 
           file_path: str, query_path: str, distance_fun: int, optimize=True) -> (dict, float):
    
    command = [
        BIN_EXE if optimize else BIN_EXE_NO_OPTIM,
        file_path,
        query_path,
        str(subseq_len), str(warp_window), str(closest_series_num), str(distance_fun)]
    
    t = time()
    stdout = subprocess.check_output(command).decode()
    t = time() - t
    parse_regexp = OUPUT_BASE_REGEXP + closest_series_num * OUPUT_STARTS_REGEXP

    closest = {}
    m_result = re.match(parse_regexp, stdout)
    if m_result:
        location, distance = m_result.groups()[:2]
        location = int(location)
        distance = float(distance)
        for loc, dist in zip(*[iter(m_result.groups()[8:])] * 2):
            closest[int(loc)] = float(dist)
    else:
        print(stdout)

    return closest, t


def search_ed(closest_series_num: int, subseq_len: int, warp_window: float, 
           data_path: str, query_path: str, normalize=True) -> (dict, float):

    """Search the most similar subseries

        The distance beetween series define as 
         ed(\{a_1, \dots, a_n\}, \{b_1, \dots, b_n\}) = \sqrt{\sum \| a_i - b_i \|_2}
    """
    
    last_location = -subseq_len
    last_dist = MAX_DIST
    closest = SortedDict()
    processing = lambda x: zscore(x, 0, ddof=0) if normalize else lambda x: x

    data = processing(np.genfromtxt(data_path))
    query = processing(np.genfromtxt(query_path)[: subseq_len])
    t = time()

    for location in range(data.shape[0] - subseq_len):
        dist = np.linalg.norm(
            processing(data[location : location + subseq_len]) - query)

        if last_location + subseq_len > location:
            if last_dist > dist:
                del closest[last_dist]
            else:
                continue

        elif closest and closest.peekitem()[0] > dist and len(closest) < closest_series_num:
            _ = closest.popitem()

        if len(closest) < closest_series_num:
            closest[dist] = location
            last_dist = dist
            last_location = location

    return dict(closest), time() - t


def decision(real_starts: list, choosen_starts: list, real_len: int,
             expected_len: int, overlap=0.8, found_twice=False) -> float:
             
    """Ratio of search subseqence
    
    Args:
        real_starts: starts of current subseries type
        choosen_starts: predicted starts
        real_len: real len of subseqences in series 
        expected_len: length of sequence that was looked for
        overlap: overlap for mark sequence as found
        found_twice: consider two choosen in one real subseq as two found
    
    Returns:
        float: ratio of the number found subsequence to all
    """

    real_starts = sorted(real_starts)
    choosen_starts = sorted(choosen_starts, reverse=True)
    real_num = len(real_starts)
    found_num = 0

    min_chosen_start = choosen_starts.pop()
    for start in real_starts:
        while True:
            if min_chosen_start is None:
                break
            if min_chosen_start < start - overlap * expected_len:
                min_chosen_start = choosen_starts.pop() if choosen_starts else None
            elif min_chosen_start > start + overlap * real_len:
                break
            else:
                found_num += 1
                min_chosen_start = choosen_starts.pop() if choosen_starts else None
                if not found_twice:
                    break

        if min_chosen_start is None:
            break

    return found_num / real_num

