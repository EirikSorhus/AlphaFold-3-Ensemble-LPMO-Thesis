from __future__ import annotations

import logging

from lpmo_pipeline.utils.logging import StructuredLogger


def test_structured_logger_reuses_single_console_handler(tmp_path) -> None:
    logger_name = "structured_logger_test"
    logger = logging.getLogger(logger_name)
    original_handlers = list(logger.handlers)
    original_propagate = logger.propagate
    for handler in list(logger.handlers):
        logger.removeHandler(handler)

    first = StructuredLogger(logger_name, tmp_path / "first")
    second = StructuredLogger(logger_name, tmp_path / "second")
    try:
        console_handlers = [
            handler
            for handler in logging.getLogger(logger_name).handlers
            if getattr(handler, StructuredLogger._CONSOLE_HANDLER_ATTR, False)
        ]
        assert len(console_handlers) == 1
        assert logging.getLogger(logger_name).propagate is False
    finally:
        first.close()
        second.close()
        for handler in list(logger.handlers):
            logger.removeHandler(handler)
        for handler in original_handlers:
            logger.addHandler(handler)
        logger.propagate = original_propagate