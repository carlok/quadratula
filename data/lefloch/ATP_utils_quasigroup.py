import itertools
import subprocess
import re

PIPE = subprocess.PIPE

MACE4 = "./mace4"
PROVER9 = "./prover9"
MAX_SECONDS = 1

### Helpers about processes
def set_max_seconds(max_seconds):
    global MAX_SECONDS
    MAX_SECONDS = max_seconds

def run_one(*popenargs, input=None, timeout=None, check=False, **kwargs):
    # todo: why not just subprocess.run?
    with subprocess.Popen(*popenargs, text=True, stdin=PIPE, stdout=PIPE, stderr=PIPE, **kwargs) as process:
        try:
            stdout, stderr = process.communicate(input, timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            process.kill()
            process.wait()
            raise
        except:
            process.kill()
            raise
        retcode = process.poll()
        if check and retcode:
            raise subprocess.CalledProcessError(retcode, process.args,
                                                output=stdout, stderr=stderr)
    return subprocess.CompletedProcess(process.args, retcode, stdout, stderr)

### Constructing tasks
def _models_task(extra0, extra1, extra2, max_models, max_seconds):
    return f'''
assign(selection_order, 2).
assign(selection_measure, 1).
assign(max_models, {max_models}).
assign(max_seconds, {max_seconds}).
assign(iterate_up_to, 100).
{extra0}
formulas(sos).
{extra1}
end_of_list.
formulas(goals).
{extra2}
end_of_list.'''

def _prove_task(extra0, extra1, extra2, max_seconds):
    return f'''
assign(sos_limit, 2000).
assign(max_weight, 1000).
assign(max_seconds, {max_seconds}).
{extra0}
formulas(sos).
{extra1}
end_of_list.
formulas(goals).
{extra2}
end_of_list.'''

### models, prove, models_or_prove
def models(*, max_models=1, max_seconds=0, check_timeout=False, extra0='', extra1='', extra2='', debug=False, timeout=None):
    if max_seconds == 0:
        max_seconds = MAX_SECONDS
    task = _models_task(extra0, extra1, extra2, max_models, max_seconds)
    if debug:
        print(task)
    try:
        maceproc = run_one(MACE4, input=task, timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        if check_timeout:
            return False
        return False
    if check_timeout and (maceproc.returncode in [1, 4, 5, 7]):
        return False
    interpretations = '\n'.join(re.findall('(?s)interpretation.*?\\]\\)\\.', maceproc.stdout))
    return [[[int(s) for s in row.split(",")]
             for row in re.split(",\n\\s*", model)]
            for model in re.findall("\\d+(?:,\\s?\\d+)+(?:,\n\\s*\\d+(?:,\\s?\\d+)+)+", interpretations)]

def prove(*, max_seconds=0, check_timeout=False, extra0='', extra1='', extra2='', timeout=None, full_output=False):
    if max_seconds == 0:
        max_seconds = MAX_SECONDS
    task = _prove_task(extra0, extra1, extra2, max_seconds)
    try:
        proverproc = run_one(PROVER9, input=task, timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        if check_timeout:
            return False
        return False
    if proverproc.returncode == 0:
        return proverproc.stdout if full_output else True
    if check_timeout and proverproc.returncode == 4:
        return False
    return False
