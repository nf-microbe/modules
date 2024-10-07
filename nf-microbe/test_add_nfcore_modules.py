
import os
import shutil

import add_nfcore_modules as add_nfcore_modules


# test changing from modules to pipeline repo
def test_pipeline_nfcore_yml():
    shutil.copy('.nf-core.yml', '.nf-core.yml.test')
    add_nfcore_modules.pipeline_nfcore_yml(nfcore_yml='.nf-core.yml.test')
    with open('.nf-core.yml.test') as f:
        content = f.read()
    assert 'repository_type: pipeline' in content
    assert 'repository_type: modules' not in content
    os.remove('.nf-core.yml.test')

# test changing from pipeline to modules repo
def test_modules_nfcore_yml():
    shutil.copy('.nf-core.yml', '.nf-core.yml.test')
    add_nfcore_modules.pipeline_nfcore_yml(nfcore_yml='.nf-core.yml.test')
    add_nfcore_modules.modules_nfcore_yml(nfcore_yml='.nf-core.yml.test')
    with open('.nf-core.yml.test') as f:
        content = f.read()
    assert 'repository_type: modules' in content
    assert 'repository_type: pipeline' not in content
    os.remove('.nf-core.yml.test')

# test installing a module
def test_install_nfcore_module():
    add_nfcore_modules.pipeline_nfcore_yml(nfcore_yml='nf-microbe/test_repo/.nf-core.yml')
    add_nfcore_modules.install_nfcore_module(directory="nf-microbe/test_repo", module_name="fastqc")
    add_nfcore_modules.modules_nfcore_yml(nfcore_yml='nf-microbe/test_repo/.nf-core.yml')
    assert os.path.exists('nf-microbe/test_repo/modules/nf-core/fastqc')
    shutil.rmtree('nf-microbe/test_repo/modules/nf-core/fastqc')
    os.remove('nf-microbe/test_repo/modules.json')

# test CLI module install
def test_main():
    add_nfcore_modules.main(["-d", "nf-microbe/test_repo", "-m", "fastqc"])
    assert os.path.exists('nf-microbe/test_repo/modules/nf-core/fastqc')
    with open('nf-microbe/test_repo/.nf-core.yml') as f:
        content = f.read()
    assert 'repository_type: modules' in content
    assert 'repository_type: pipeline' not in content
    shutil.rmtree('nf-microbe/test_repo/modules/nf-core/fastqc')
    os.remove('nf-microbe/test_repo/modules.json')
