import os
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML

# requirements
# ------------
# brew install cairo
# brew install pango


def print_pdf(df, title, filename):
    """

    :param str filename: name of pdf file to write
    :param str title: title for pdf table
    :param pd.DataFrame df: data to be printed to pdf
    """

    # load template
    module_dir = os.path.split(__file__)[0]
    env = Environment(loader=FileSystemLoader(module_dir + '/html_/'))
    template = env.get_template('pdf_table.html')

    # check that filename has correct suffix
    if not filename.endswith('.pdf'):
        filename += '.pdf'

    # generate html
    template_vars = {
        'title': title,
        'data_frame': df.to_html()
    }
    html_out = template.render(template_vars)

    # create css styled pdf table from html
    HTML(string=html_out).write_pdf(filename, stylesheets=[module_dir + '/css/style.css'])
