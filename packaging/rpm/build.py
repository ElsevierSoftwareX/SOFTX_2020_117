import os
import textwrap

from subprocess import call, check_output

from collections import defaultdict
import string

class SafeDict(dict):
    def __missing__(self, key):
        return '{' + key + '}'

os.chdir("/root/pykat")

call("git pull".split())

git_describe = str(check_output('git describe --tags'.split())).split("-")

os.chdir("/root/pykat")

vals = {
    "version": git_describe[0],
    "release": git_describe[1],
}

with open("/root/pykat/packaging/rpm/pykat.spec", "w") as f:
    s = string.Formatter().vformat(textwrap.dedent("""
            %define name pykat
            %define version {version}
            %define unmangled_version {version}.{release}
            %define release {release}

            Summary: Python interface and tools for FINESSE
            Name: %{name}
            Version: %{version}
            Release: %{release}
            Source0: %{name}-%{unmangled_version}.tar.gz
            License: GPL v2
            Group: Development/Libraries
            BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-buildroot
            Prefix: %{_prefix}
            BuildArch: noarch
            Vendor: Daniel Brown <finesse@star.sr.bham.ac.uk>
            Url: http://gwoptics.org/pykat

            %description

            %prep
            %setup -n %{name}-%{unmangled_version} -n %{name}-%{unmangled_version}

            %build
            python3 setup.py build

            %install
            python3 setup.py install --single-version-externally-managed -O1 --root=$RPM_BUILD_ROOT --record=INSTALLED_FILES

            %clean
            rm -rf $RPM_BUILD_ROOT

            %files -f INSTALLED_FILES
            %defattr(-,root,root)
            """), (), SafeDict(**vals))
            
    f.write(s)

call("python setup.py bdist_rpm".split())

#call("cp /root/pykat/packaging/rpm/rpmbuild/RPMS/x86_64/finesse-{version}-{release}.x86_64.rpm /host".format(**vals).split())


