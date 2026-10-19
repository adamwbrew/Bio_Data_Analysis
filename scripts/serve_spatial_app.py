#!/usr/bin/env python3
"""Open the saved spatial explorer in a browser, without serving the repository."""
import argparse
import errno
import json
from pathlib import Path
import sys
from urllib.error import HTTPError, URLError
from urllib.request import urlopen
import webbrowser

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from bio_visuals.standalone_app import DEFAULT_PORT, HEALTH_PATH, create_server, document_identity


def port_number(value):
    number = int(value)
    if not 1 <= number <= 65535:
        raise argparse.ArgumentTypeError("Port must be between 1 and 65535.")
    return number


def same_app_running(url, app_file):
    try:
        with urlopen(url.rstrip("/") + HEALTH_PATH, timeout=2) as response:
            identity = json.loads(response.read(2048))
        return identity == document_identity(app_file)
    except (HTTPError, URLError, TimeoutError, OSError, ValueError):
        return False


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=port_number, default=DEFAULT_PORT,
                        help="Local port (default: %(default)s).")
    parser.add_argument("--no-browser", action="store_true", help="Print the URL without opening a browser.")
    args = parser.parse_args(argv)
    app_file = ROOT / "outputs" / "spatial_handoff" / "app" / "index.html"
    url = "http://127.0.0.1:{}/".format(args.port)
    try:
        server = create_server(app_file, args.port)
    except FileNotFoundError as error:
        parser.exit(1, str(error) + "\n")
    except OSError as error:
        if error.errno == errno.EADDRINUSE:
            if same_app_running(url, app_file):
                print("The same spatial explorer is already running: " + url, flush=True)
                if not args.no_browser:
                    webbrowser.open(url)
                return 0
            parser.exit(1, "Port {} is already used by another app or snapshot. "
                        "Choose another port, for example --port {}.\n".format(
                            args.port, args.port + 1 if args.port < 65535 else 8050))
        parser.exit(1, "Could not start the local spatial explorer: {}\n".format(error))
    print("PLA–ChAT Spatial Explorer: " + url, flush=True)
    print("Local computer only. Keep this window open; Control+C stops the app.", flush=True)
    if not args.no_browser:
        webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nSpatial explorer stopped.")
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
