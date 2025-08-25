import os
import pytest
from kimmdy import cmd


def test_get_cmdline_args(mocker):
    parser = mocker.MagicMock(spec=cmd.argparse.ArgumentParser())
    mocker.patch("kimmdy.cmd.argparse.ArgumentParser", autospec=True).return_value = (
        parser
    )
    cmd.get_cmdline_args()
    parser.add_argument.assert_called()


def test_run_plugins(mocker):
    mock_args = mocker.sentinel
    mock_args.show_plugins = True
    with pytest.raises(SystemExit):
        cmd._run(mock_args)


def test_run_no_input(mocker):
    mock_args = mocker.sentinel
    mock_args.show_plugins = False
    mock_args.input = "./DoesNotExists"
    with pytest.raises(FileNotFoundError):
        cmd._run(mock_args)


def test_run_write_jobscript(mocker, tmp_path):
    mock_args = mocker.sentinel
    mock_args.show_plugins = False
    mock_args.input = "./"
    mock_args.generate_jobscript = True

    mocker.patch("kimmdy.cmd.Config")
    mocker.patch("kimmdy.cmd.RunManager")

    os.chdir(tmp_path)
    cmd._run(mock_args)
    assert len(list(tmp_path.glob("jobscript*"))) == 1


def test_entry_point(mocker):
    ga = mocker.patch("kimmdy.cmd.get_cmdline_args")
    run = mocker.patch("kimmdy.cmd._run")
    cmd.entry_point_kimmdy()
    ga.assert_called_once()
    run.assert_called_once()
