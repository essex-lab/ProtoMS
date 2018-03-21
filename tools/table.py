"""Collection of classes to handle output table formatting"""

from itertools import izip_longest


class Table(object):
    """Output table for displaying data"""
    def __init__(self, title, fmts, headers=[], space=3):
        """Parameters
        ----------
        title: string
            title printed above the output table
        fmts: list of strings
            valid format strings e.g. '%.3f' for each column's entries
        headers: list of strings
            top entry of each column
        space: int
            number of spaces between columns
        """
        self.title = title
        self.columns = [
            Column(head, fmt, fmt)
            for fmt, head in izip_longest(fmts, headers, fillvalue='')
        ]
        self.spacer = " " * space

    def add_row(self, data):
        """Add a row to the table.

        Parameters
        ----------
        data: list of objects
            table entries for this run, objects must be convertable to
            strings according to the corresponding column format
        """
        for col, dat in zip(self.columns, data):
            col.add_data(dat)

    def add_blank_row(self):
        """Add a blank row to the table."""
        for col in self.columns:
            col.add_data(None)

    def __str__(self):
        s = ''
        if self.title:
            s += "%s\n" % self.title

        for col in self.columns:
            s += ("%-{}s" + self.spacer).format(col.width) % col.header
        s += '\n'

        for row in zip(*self.columns):
            s += self.spacer.join(row) + '\n'
        return s


class Column(object):
    """Table column entry."""
    def __init__(self, header='', value_fmt="%.4f", error_fmt="%.4f"):
        """Parameters
        ----------
        header: string
            top row descriptive entry
        value_fmt: string
            format string for the value component of each entry
        error_fmt: string
            format string for the error component of each entry
        """
        self.header = header
        self.left = SubColumn(value_fmt)
        self.right = SubColumn(error_fmt)

    def add_data(self, data):
        """Add a entry to this column.

        Parameters
        ----------
        data: Quantity object or number
            entry to add
        """
        try:
            self.left.add_data(data.value)
            self.right.add_data(data.error)
        except AttributeError:
            self.left.add_data(data)
            self.right.add_data(None)

    def __iter__(self):
        for l, r in zip(self.left, self.right):
            stuff = "%s +- %s" % (l, r) if len(r) != r.count(" ") else "%s" % l
            yield "%-{}s".format(self.width) % stuff

    @property
    def width(self):
        """Width of the column in number of characters"""
        if self.right.width:
            contents = self.left.width + self.right.width + 4
        else:
            contents = self.left.width
        return max(contents, len(self.header))


class SubColumn(object):
    """The left or right side of a Column."""
    def __init__(self, fmt="%.4f"):
        """Parameters
        ----------
        fmt: string, optional
            format string for the entries of this Subcolumn
        """
        self.values = []
        self.fmt = fmt
        self.width = 0

    def add_data(self, data):
        """Add a entry to this subcolumn.

        Parameters
        ----------
        data: number
            entry to add
        """
        line = self.fmt % data if data is not None else ''
        self.values.append(line)
        self.width = max(len(line), self.width)

    def __iter__(self):
        for val in self.values:
            yield "%{}s".format(self.width) % val
