"""
Functions to multithread process a list of inputs.
"""

from queue import Queue
import threading

from pdm_utils.classes.progressbar import ProgressBar, show_progress


class MixedThread(threading.Thread):
    def __init__(self, thread_id, work_queue, work_lock,
                 result_queue, result_lock):
        threading.Thread.__init__(self)

        self.name = thread_id

        # Queue filled with tuples containing target function args
        self.work_queue = work_queue
        self.result_queue = result_queue
        # Lock controlling access to queue accross multiple threads
        self.work_lock = work_lock
        self.result_lock = result_lock

    def run(self):
        while not self.work_queue.empty():
            # Tell other threads "I'm accessing the queue right now"
            self.work_lock.acquire()

            job = self.work_queue.get()
            # Let other threads have a turn
            self.work_lock.release()

            target = job[0]
            work_item = job[1]

            # Run target function with args retrieved from queue
            result = target(*work_item)

            # Repeat with the result queue/lock
            self.result_lock.acquire()
            self.result_queue.put(result)
            self.result_lock.release()


def create_threads(work_queue, work_lock, result_queue, result_lock,
                   num_threads):
    """Function to create a pseudo-MixIn Thread class for multithreading.

    :param target: Function used by the thread to process work items.
    :type target: Function
    :param work_queue: Queue containing tuples of target function args.
    :type work_queue: Queue
    :param work_lock: Lock for controlling access of work queue accross threads
    :type work_lock: Lock
    :param result_queue: Queue to place results into
    :type result_queue: Queue
    :param result_lock: Lock for controlling result queue accross threads.
    :type result_lock: Lock
    :param num_threads: Number of threads to be created to process work stack.
    :type num_threads: int
    :returns: Returns a list of MixIn Thread objects with target function
    :rtype: list[Thread]
    """
    threads = []
    for x in range(num_threads):
        threads.append(MixedThread(x+1, work_queue, work_lock,
                                   result_queue, result_lock))

    return threads


def multithread(work_items, threads, target, verbose=False):
    """Runs list of work items with specified number of threads.

    :param target: Function used by the thread to process work items.
    :type target: Function
    :type work_items: List containing tuples of target function args.
    :type work_stack: list
    :param threads: Number of threads to be created to process work_stack.
    :type threads: int
    """
    # Create lock which prevents work queue access deadlocks
    work_lock = threading.Lock()
    work_queue = Queue()

    # Create lock which prevents result queue access deadlocks
    result_lock = threading.Lock()
    result_queue = Queue()

    tasks = 0
    if verbose:
        interval = max([1, len(work_items) // 100])

        for i in range(len(work_items)):
            if i % interval == 0:
                tasks += 1
                work_queue.put((show_progress, (i, len(work_items))))

            tasks += 1
            work_queue.put((target, work_items[i]))

        tasks += 1
        work_queue.put((show_progress, (len(work_items), len(work_items))))
    else:
        for item in work_items:
            tasks += 1
            work_queue.put((target, item))

    threads = create_threads(work_queue, work_lock, result_queue, result_lock,
                             threads)

    for thread in threads:
        # Calls thread run()
        thread.start()

    for thread in threads:
        # Blocks and waits for threads to finish
        thread.join()

    # Remove non-Progress results
    results = []
    for i in range(tasks):
        result = result_queue.get()
        if result is not None and not isinstance(result, ProgressBar):
            results.append(result)

    return results
