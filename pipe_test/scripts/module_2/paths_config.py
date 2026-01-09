from pathlib import Path
from datetime import datetime


class PipelinePaths:
    base_dir = Path(".")
    metadata_dir = base_dir / "data" / "metadata"
    run_dir = base_dir / "data" / "run"
    sequences_dir = base_dir / "data" / "sequences"
    logs_dir = base_dir / "logs"
    cazy_dir = base_dir / "cazy_raw"

    @classmethod
    def configure_base(cls, base_dir: Path) -> None:
        cls.base_dir = Path(base_dir).resolve()
        cls.metadata_dir = cls.base_dir / "data" / "metadata"
        cls.run_dir = cls.base_dir / "data" / "run"
        cls.sequences_dir = cls.base_dir / "data" / "sequences"
        cls.logs_dir = cls.base_dir / "logs"
        cls.cazy_dir = cls.base_dir / "cazy_raw"

    @staticmethod
    def timestamp():
        return datetime.now().strftime("%Y%m%d_%H%M%S")

    @classmethod
    def make_all_dirs(cls):
        for d in [cls.metadata_dir, cls.run_dir, cls.sequences_dir, cls.logs_dir, cls.cazy_dir]:
            d.mkdir(parents=True, exist_ok=True)
