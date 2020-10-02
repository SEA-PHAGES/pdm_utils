"""
Functions to multithread process a list of inputs.
"""

from queue import LifoQueue
import threading


def create_threads(target, work_stack, stack_lock, num_threads):
    """Function to create a pseudo-MixIn Thread class for multithreading.

    :param target: Function used by the thread to process work items.
    :type target: Function
    :param work_stack: Stack containing tuples of target function args.
    :type work_stack: LifoQueue
    :param stack_lock: Lock for controlling access of stack accross threads.
    :type stack_lock: Lock
    :param num_threads: Number of threads to be created to process work stack.
    :type num_threads: int
    :returns: Returns a list of MixIn Thread objects with target function
    :rtype: list[Thread]
    """
    class MixedThread(threading.Thread):
        def __init__(self, threadid):
            threading.Thread.__init__(self)

            self.name = threadid

            self.target = target
            self.stack = work_stack
            self.lock = stack_lock

        def run(self):
            while not self.stack.empty():
                self.lock.acquire()

                work_item = self.stack.get()
                self.lock.release()

                self.target(*work_item)

    threads = []
    for x in range(num_threads):
        threads.append(MixedThread(x+1))

    return threads


def multithread(target, work_items, threads=1):
    """Runs list of work items with specified number of threads.

    :param target: Function used by the thread to process work items.
    :type target: Function
    :type work_items: List containing tuples of target function args.
    :type work_stack: list
    :param threads: Number of threads to be created to process work_stack.
    :type threads: int
    """
    lock = threading.Lock()
    stack = LifoQueue()

    for item in work_items:
        stack.put(item)

    threads = create_threads(target, stack, lock, threads)

    for thread in threads:
        thread.start()

    for thread in threads:
        thread.join()
