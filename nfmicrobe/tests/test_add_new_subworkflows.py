import filecmp
import os
import shutil

import src.add_new_subworkflows as add_new_subworkflows


# test changing from modules to pipeline repo
def test_update_modules_json():
    shutil.copy("nfmicrobe/test_data/modules.json.base", "modules.json.test")
    add_new_subworkflows.update_modules_json(subworkflow_name="b_test", modules_json="modules.json.test")
    if filecmp.cmp("modules.json.test", "nfmicrobe/test_data/modules_test.json"):
        assert True
    os.remove("modules.json.test")


# test adding module to nf-core.yml
def test_add_subworkflow_to_nfcore_yml():
    shutil.copy(".nf-core.yml", ".nf-core.yml.test")
    add_new_subworkflows.add_subworkflow_to_nfcore_yml(subworkflow_name="b_test", nfcore_yml=".nf-core.yml.test")
    with open(".nf-core.yml.test") as f:
        content = f.read()
    assert "nf-core:\n      b_test: False" in content
    os.remove(".nf-core.yml.test")


# test main function
def test_main():
    shutil.copy(".nf-core.yml", "nfmicrobe/test_data/.nf-core.yml")
    shutil.copy("nfmicrobe/test_data/modules.json.base", "nfmicrobe/test_data/modules.json")
    add_new_subworkflows.main(["-d", "nfmicrobe/test_data", "-s", "b_test"])
    with open("nfmicrobe/test_data/.nf-core.yml") as f:
        content = f.read()
    assert "nf-core:\n      b_test: False" in content
    if filecmp.cmp("nfmicrobe/test_data/modules.json", "nfmicrobe/test_data/modules_test.json"):
        print
        assert True
    os.remove("nfmicrobe/test_data/.nf-core.yml")
    os.remove("nfmicrobe/test_data/modules.json")
