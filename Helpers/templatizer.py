import os
import string
import argparse

__author__ = "Carles Perez Lopez"


def parse_arguments():
    """
        Parse user arguments
        Output: list with all the user arguments
    """

    parser = argparse.ArgumentParser(description="""Description: Program to replace a pattern from a templatized file
    and write several new files with this pattern replaced by the different elements of a list.""")
    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument("template", type=str, help="Template file.")
    required_named.add_argument("--subs_list", nargs="+",
                                help="List of elements to add and replace by the pattern. They will be appended.")
    parser.add_argument("--pattern", type=str, default="NAME",
                        help="Pattern to be replaced. Default: 'NAME'. In template: $NAME")
    args = parser.parse_args()

    return args.template, args.subs_list, args.pattern


def read_file(file_path):
    """
    Reads an input file and returns an string
    :param file_path: path to the input file
    :return: content of the file
    """
    with open(file_path, "r") as f:
        content = f.read()
    return content


def write_file(content, out_file_path):
    """
    Writes a string to an output file.
    :param content: content of the output file
    :type content: str
    :param out_file_path: output file path.
    :return: content of the file
    """
    with open(out_file_path, "w") as f:
        f.write(content)


def add_suffix_to_file(file_path, suffix):
    """
    Given a filename it adds a suffix before the extension.
    :param file_path:
    :param suffix:
    :return:
    """
    name, ext = os.path.splitext(file_path)
    name = name + suffix
    return "".join([name, ext])


def define_subs(pattern, value):
    key_dict = {pattern: value}
    return key_dict


def replace_list(template, pattern, list_of_replacements):
    """
    Given a template file, it replaces the selected pattern by the different elements of the list of replacements. Each
    element of the list writes a new output file with its new content replaced.
    :param template: templatized input file
    :param pattern: pattern to be replaced (example: $SOMETHING)
    :param list_of_replacements: list of elements to replace for the pattern
    :return: writes as many files as elements in the list of replacements, using also their names as suffix for
    output files.
    """
    temp_content = read_file(template)
    for sub in list_of_replacements:
        t = string.Template(temp_content)
        keys = define_subs(pattern=pattern, value=sub)
        replaced_content = t.substitute(keys)
        output = add_suffix_to_file(template, sub)
        write_file(replaced_content, output)


if __name__ == '__main__':
    template, subs_list, pattern = parse_arguments()
    replace_list(template, pattern, subs_list)


