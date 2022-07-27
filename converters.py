"""Convert between data types."""
import argparse
import configparser

def argparser_to_configparser(
    args: argparse.Namespace,
    config: configparser.ConfigParser,
    section: str,
) -> configparser.ConfigParser:
    """Update a `configparser` object from the values in `args`.

    If option is in `args` but not in `config`, then it is skipped.

    Args:
      args: argument parser object.
      config: configuration parser.
      section: section of `config` to update

    Returns:
      An updated configuration parser.
    """
    # Initial values
    new_values = {}
    args_dict = vars(args)

    # Filter values
    for key, val in args_dict.items():
        if key in config.options(section):
            new_values[key] = val

    # Update config parser
    config.read_dict({section: new_values})

    return config
