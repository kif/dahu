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
    server_instance="dahu",
    domain="id00",
    family="dahu",
    member="1",
    device_class="DahuDS",
):
    """Register a TANGO device of particular device class"""
    dev_name = "/".join([domain, family, member])
    db = PyTango.Database()
    try:
        proxy = PyTango.DeviceProxy(dev_name)
    except PyTango.DevFailed:
        dev_info = PyTango.DbDevInfo()
        dev_info.name = dev_name
        dev_info._class = device_class
        dev_info.server = "/".join([device_class, server_instance])
        db.add_device(dev_info)
        proxy = PyTango.DeviceProxy(dev_name)
        logger.info("Registered a new TANGO device '%s'", get_uri(proxy))
    else:
        dev_info = db.get_device_info(dev_name)
        server_instance2 = dev_info.ds_full_name.split("/")[-1]
        if server_instance == server_instance2:
            logger.info("TANGO device '%s' is already registered", get_uri(proxy))
        else:
            logger.info(
                "TANGO device '%s' is already registered with another server '%s'",
                get_uri(proxy),
                server_instance2,
            )
            server_instance = server_instance2
    logger.info(
        "To start the TANGO device server:\n\n   TANGO_HOST=%s:%s dahu-server %s\n",
        proxy.get_db_host(),
        proxy.get_db_port(),
        server_instance,
    )
    logger.info(
        "To get a proxy to the device:\n\n   proxy = DeviceProxy('%s')\n",
        get_uri(proxy),
    )


def main(argv=None):
    if argv is None:
        argv = sys.argv

    description = """Register a Dahu TANGO device

    The full device name "<domain>/<family>/<member>" ("id00/dahu/1" by default)
    needs to be unique in the scope of one TANGO database.

    TANGO devices are instantiated in TANGO device servers ("dahu" by default).

    Each TANGO device server corresponds to one system process.
    """

    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--server",
        type=str,
        default="dahu",
        help="Server instance name (Default: 'dahu')",
    )
    parser.add_argument(
        "--domain",
        type=str,
        default="id00",
        help="Device domain name (Default: 'id00')",
    )
    parser.add_argument(
        "--family",
        type=str,
        default="dahu",
        help="Device family name (Default: 'dahu')",
    )
    parser.add_argument(
        "--member",
        type=str,
        default="1",
        help="Device name (Default: '1')",
    )
    args = parser.parse_args(argv[1:])
    register(
        server_instance=args.server,
        domain=args.domain,
        family=args.family,
        member=args.member,
    )


if __name__ == "__main__":
    sys.exit(main())
