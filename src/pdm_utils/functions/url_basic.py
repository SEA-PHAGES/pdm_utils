import re
import sys
import urllib3


# GET FUNCTIONS
# -----------------------------------------------------------------------------
def pool_request(url, pool=False,
                 pipeline=False, preload=False, expect_status=200):
    """
    Create a urllib3 PoolManager to access the url link for list of databases
    :returns An HTTPResponse object
    :rtype: urllib3.response.HTTPResponse
    """
    if pool is None:
        pool = urllib3.PoolManager()
    response = pool.request("GET", url, preload_content=preload)

    pool.clear()

    if pipeline:
        if response.status != expect_status:
            print("Received invalid response from server.\n"
                  "Aborting pipeline.")
            sys.exit(1)

    return response


def download_file(file_url, filepath, pool=None, chunk_size=512000,
                  preload=False, verbose=True):
    """
    Retrieve a file from the server.
    :param file_url:  URL for the desired file:
    :type file_url: str
    :param filepath: Local path where the file can be downloaded to.
    :type filepath: Path
    :returns: Status of the file retrieved from the server
    :rtype: bool
    """
    request = pool_request(file_url, pool=pool, pipeline=False)
    status = request.status

    if status == 200:
        filesize = int(request.getheader("Content-Length"))
        output_file = filepath.open(mode="wb")
        for chunk in request.stream(chunk_size):
            output_file.write(chunk)
            progress = int(request._fp_bytes_read / filesize * 100)
            if verbose:
                print("\r[{}{}] {}%".format("#" * int(progress / 2),
                                            " " * (50 - int(progress / 2)),
                                            progress), end="")
    else:
        if verbose:
            print(" ".join(["ERROR.  HTTP Response", str(status)]))

    if verbose:
        print("")

    return status


def get_url_listing_data(response, get_dates=False):
    """
    Get list of files and directories from the link using a response object
    :param response: PoolManager response for the specified url link
    :type response: urllib3.response.HTTPResponse
    :returns: A list of entry listings
    :rtype: list
    """
    listings = list()

    response_data = response.data.decode("utf-8")
    # print(response_data)

    # Get database names
    names_regex = """<a href="([-_\w\d./]+)">"""
    names_regex = re.compile(names_regex)
    names = names_regex.findall(response_data)

    if get_dates:
        # Get database dates
        dates_regex = ("""<td align="right">"""
                       "(\d+[-]\d+[-]\d+)\s+\d+[:]\d+\s+</td>")
        dates_regex = re.compile(dates_regex)
        dates = dates_regex.findall(response_data)

    response.close()

    # creating a list of dictionaries
    for i in range(len(names)):
        db_dict = dict()
        db_dict["num"] = i+1
        db_dict["name"] = names[i]
        if get_dates:
            db_dict["date"] = dates[i]

        listings.append(db_dict)

    return listings


# PARSING FUNCTIONIS
# ----------------------------------------------------------------------------
def get_url_listing_names(response, suffix=None):
    listing = get_url_listing_data(response)

    listing_names = list()

    for data_dict in listing:
        name = data_dict.get("name")

        if name is None:
            continue

        if suffix:
            if len(suffix) >= len(name):
                continue

            if name.endswith(suffix):
                listing_names.append(name[:-(len(suffix))])
        else:
            listing_names.append(name)

    return listing_names


def get_url_listing_dirs(response):
    return get_url_listing_names(response, suffix="/")


def get_url_listing_files(response, file_ext):
    return get_url_listing_names(response, "".join([".", str(file_ext)]))
