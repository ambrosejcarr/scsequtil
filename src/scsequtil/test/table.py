import numpy as np
import pandas as pd
import unittest
import scsequtil.table


class TestTable(unittest.TestCase):

    def test_table(self):
        """create a table from some mock data

        :return:
        """
        data = np.zeros(4).reshape(2, 2)
        df = pd.DataFrame(data, index=['green', 'eggs'], columns=['foo', 'bar'])
        scsequtil.table.print_pdf(df, 'a table', 'table.pdf')


if __name__ == "__main__":
    unittest.main()
