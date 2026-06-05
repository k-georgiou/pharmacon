"""
Pharmacon: A Molecular Dynamics Simulation Analysis Toolkit
    Copyright© 2026  Kyriakos Georgiou

Tests for pharmacon.command_line.registry.CommandRegistry
"""
import pytest

from pharmacon.command_line.registry import CommandRegistry, get_registry
from pharmacon.command_line.base import CommandSpec, SubcommandSpec



class TestCommandRegistryDiscovery:
    @pytest.fixture(scope="module")
    def registry(self):
        return CommandRegistry()

    def test_registry_is_non_empty(self, registry):
        assert len(registry.commands) > 0

    @pytest.mark.parametrize("name", ["trajectory", "structure", "merge", "plot", "export", "dump"])
    def test_expected_command_exists(self, registry, name):
        assert name in registry.commands, f"Expected command '{name}' not discovered"

    def test_get_command_returns_command_spec(self, registry):
        cmd = registry.get_command("trajectory")
        assert isinstance(cmd, CommandSpec)

    def test_get_command_has_name_field(self, registry):
        cmd = registry.get_command("trajectory")
        assert cmd.name == "trajectory"

    def test_get_command_unknown_returns_none(self, registry):
        assert registry.get_command("absolutely-nonexistent") is None

    def test_commands_property_is_dict(self, registry):
        assert isinstance(registry.commands, dict)


class TestSubcommandDiscovery:
    @pytest.fixture(scope="module")
    def registry(self):
        return CommandRegistry()

    @pytest.mark.parametrize("subcommand", ["rmsd", "distances", "angles"])
    def test_trajectory_subcommands_exist(self, registry, subcommand):
        sub = registry.get_subcommand("trajectory", subcommand)
        assert sub is not None, f"Expected trajectory subcommand '{subcommand}'"

    def test_subcommand_is_subcommand_spec(self, registry):
        sub = registry.get_subcommand("trajectory", "rmsd")
        assert isinstance(sub, SubcommandSpec)

    def test_subcommand_has_callable_run_fn(self, registry):
        sub = registry.get_subcommand("trajectory", "rmsd")
        assert callable(sub.run_fn)

    def test_subcommand_has_callable_build_parser_fn(self, registry):
        sub = registry.get_subcommand("trajectory", "rmsd")
        assert callable(sub.build_parser_fn)

    def test_subcommand_has_callable_validate_fn(self, registry):
        sub = registry.get_subcommand("trajectory", "rmsd")
        assert callable(sub.validate_fn)

    def test_subcommand_has_non_empty_summary(self, registry):
        sub = registry.get_subcommand("trajectory", "rmsd")
        assert sub.summary

    def test_get_subcommand_unknown_command_returns_none(self, registry):
        assert registry.get_subcommand("nonexistent", "rmsd") is None

    def test_get_subcommand_unknown_subcommand_returns_none(self, registry):
        assert registry.get_subcommand("trajectory", "nonexistent") is None


class TestGetRegistrySingleton:
    def test_get_registry_returns_command_registry(self):
        r = get_registry()
        assert isinstance(r, CommandRegistry)

    def test_get_registry_returns_same_instance(self):
        r1 = get_registry()
        r2 = get_registry()
        assert r1 is r2

    def test_singleton_has_commands(self):
        r = get_registry()
        assert len(r.commands) > 0
