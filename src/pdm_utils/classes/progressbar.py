class ProgressBar:
    def __init__(self, current, end, width=50):
        self._percent = int(float(current) / end * 100)
        self._width = int(width)
        self._width_ratio = int(100 / width)
        self._multiplier = int(self._percent / self._width_ratio)
        self._padding = self._width - self._multiplier

    def __str__(self):
        s = f"[{'#' * self._multiplier}{' ' * self._padding}] {self._percent}%"
        return s


def show_progress(current, end, width=50):
    """
    Creates a ProgressBar instance and prints it in-line if
    current < end. If current == end, prints with newline.
    :param current: current step (1 through n)
    :type current: int
    :param end: number of steps (n)
    :type end: int
    :param width: character width for the progressbar
    :type width: int
    :return: progressbar
    """
    progressbar = ProgressBar(current, end, width)
    if current < end:
        print(f"\r{progressbar}", end="")
    else:
        print(f"\r{progressbar}")
