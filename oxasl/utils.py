"""
Misc utility functions
"""

import six

class Tee(object):
    """
    Output stream which keeps a string record of everything
    written to it, and also sends it on to any number of
    sub-streams
    """

    def __init__(self, *streams):
        self._streams = [six.StringIO(),]
        self._streams.extend(streams)

    def add(self, stream):
        """
        Add a sub-stream
        """
        if stream:
            self._streams.append(stream)

    def write(self, text):
        """
        Write to all output streams
        """
        for stream in self._streams:
            stream.write(text)

    def flush(self):
        """
        Flush all output streams
        """
        for stream in self._streams:
            stream.flush()

    def __str__(self):
        return self._streams[0].getvalue()
