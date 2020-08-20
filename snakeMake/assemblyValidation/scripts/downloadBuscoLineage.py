from busco.BuscoConfig import BuscoConfigMain
from busco.ConfigManager import BuscoConfigManager
from busco.BuscoLogger import BuscoLogger
import sys

usage = "python3 {} <input lineage_name>".format(sys.argv[0])

if len(sys.argv) != 2:
    print(usage)
    sys.exit()

params = {'config_file' : None}

manager = BuscoConfigManager(params)

configfile = manager.get_config_file()

cobj = BuscoConfigMain(configfile, params, sys.argv)

lineage_dataset = cobj.get("busco_run", "lineage_dataset")
cobj.set_results_dirname(lineage_dataset)
cobj.config.download_lineage_file(lineage_dataset)
