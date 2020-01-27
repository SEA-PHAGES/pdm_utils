"""
Functions to parallelize of processing of a list of inputs.
Adapted from https://docs.python.org/3/library/multiprocessing.html
"""

import multiprocessing as mp

from pdm_utils.functions.basic import show_progress


def parallelize(inputs, num_processors, task):
    """
    Parallelizes some task on an input list across the specified number
    of processors
    :param inputs: list of inputs
    :param num_processors: number of processor cores to use
    :param task: name of the function to run
    :return: results
    """
    results = []

    # Don't do any work if there are no inputs
    if len(inputs) == 0:
        return results

    num_processors = count_processors(inputs, num_processors)

    tasks = []
    for item in inputs:
        if not isinstance(item, tuple):
            item = (item,)
        tasks.append((task, item))

    # Start working on the jobs
    results = start_processes(tasks, num_processors)

    return results


def count_processors(inputs, num_processors):
    """
    Programmatically determines whether the specified num_processors is
    appropriate. There's no need to use more processors than there are
    inputs, and it's impossible to use fewer than 1 processor or more
    than exist on the machine running the code.
    :param inputs: list of inputs
    :param num_processors: specified number of processors
    :return: num_processors (optimized)
    """
    if num_processors < 1 or num_processors > mp.cpu_count():
        print(f"Invalid number of CPUs specified ({num_processors})")
        num_processors = min([mp.cpu_count(), len(inputs)])
        print(f"Using {num_processors} CPUs...")

    return num_processors


def worker(input_queue, output_queue):
    for func, args in iter(input_queue.get, 'STOP'):
        result = func(*args)
        output_queue.put(result)
    return


def start_processes(inputs, num_processors):
    """
    Creates input and output queues, and runs the jobs
    :param inputs: jobs to run
    :param num_processors: optimized number of processors
    :return: results
    """
    # Make sure new processes are forked, not spawned
    mp.set_start_method("fork")

    job_queue = mp.Queue()
    done_queue = mp.Queue()

    # Counter so we know how many tasks we have in all (input + show_progress tasks)
    task_count = 0

    # Put inputs into job queue
    for i in range(len(inputs)):
        task_count += 1
        if i % 10 == 0:
            job_queue.put((show_progress, (i, len(inputs))))
            task_count += 1
        job_queue.put(inputs[i])
    task_count += 1
    job_queue.put((show_progress, (len(inputs), len(inputs))))

    # Put a bunch of 'STOP' signals at the end of the queue
    for i in range(num_processors):
        job_queue.put('STOP')

    # Start up workers
    worker_pool = []
    for i in range(num_processors):
        worker_n = mp.Process(target=worker, args=(job_queue, done_queue))
        worker_n.start()
        worker_pool.append(worker_n)

    # Grab results from done queue
    results = []

    # Remove non-list results (e.g. progress results)
    for i in range(task_count):
        result = done_queue.get()
        if isinstance(result, list):
            results.append(result)

    [worker_n.join() for worker_n in worker_pool]

    # Leave the progress bar line
    print("\n")

    return results
