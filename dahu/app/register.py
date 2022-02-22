import sys
import logging
import argparse
import PyTango

logger = logging.getLogger("dahu_server")


def get_uri(proxy):
    return "tango://{}:{}/{}".format(
        proxy.get_db_host(), proxy.get_db_port(), proxy.dev_name()
    )


logging.basicConfig(level=logging.INFO)


def register(
    server="DahuDS", instance="dahu", domain="id00", family="dahu", member="1"
):
    """Register TANGO device of class DahuDSClass"""
    dev_name = "/".join([domain, family, member])
    try:
        proxy = PyTango.DeviceProxy(dev_name)
        logger.info("'%s' already registered", get_uri(proxy))
    except PyTango.DevFailed:
        db = PyTango.Database()
        dev_info = PyTango.DbDevInfo()
        dev_info.name = dev_name
        dev_info._class = "DahuDSClass"
        server = "/".join([server, instance])
        dev_info.server = server
        db.add_device(dev_info)
        proxy = PyTango.DeviceProxy(dev_name)
        logger.info("'%s' registered", get_uri(proxy))
    logger.info(
        "To start the server:\n\n   TANGO_HOST=%s:%s dahu-server %s\n",
        proxy.get_db_host(),
        proxy.get_db_port(),
        instance,
    )
    logger.info("To get a proxy to the device:\n\n   DeviceProxy('%s')\n", get_uri(proxy))


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(description="Register a Dahu Tango device")
    parser.add_argument(
        "--instance",
        type=str,
        default="dahu",
        help="Server instance name (session name by default)",
    )
    parser.add_argument(
        "--domain",
        type=str,
        default="id00",
        help="Device domain name",
    )
    parser.add_argument(
        "--family",
        type=str,
        default="dahu",
        help="Device family name",
    )
    parser.add_argument("--member", type=str, default="1", help="Device name")
    args = parser.parse_args(argv[1:])
    register(
        instance=args.instance,
        domain=args.domain,
        family=args.family,
        member=args.member,
    )


if __name__ == "__main__":
    sys.exit(main())
